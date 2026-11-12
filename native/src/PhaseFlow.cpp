#include "PhaseSpaceInternal.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <limits>
namespace dynsyslab::core::phase_detail {
double norm(const PhaseState &x, const PhaseSpaceOptions &o) {
    double r = 0;
    for (std::size_t j = 0; j < x.size(); ++j)
        r = std::max(r, std::abs(x[j]) / (o.upper[j] - o.lower[j]));
    return r;
}
bool inside(const PhaseState &x, const PhaseSpaceOptions &o) {
    for (std::size_t j = 0; j < o.lower.size(); ++j)
        if (!std::isfinite(x[j]) || x[j] < o.lower[j] || x[j] > o.upper[j])
            return false;
    return true;
}
PhaseState fieldAt(const ODESystem &system, const PhaseState &x, double time, bool periodic) {
    const auto n = system.descriptor().stateDimension;
    PhaseState f(n);
    if (x.size() < n || system.invalid(time, x.data()))
        throw std::domain_error("State is outside the system domain");
    if (periodic) system.rhs(time, x.data(), f.data());
    else system.rhsWithInput(0, x.data(), 0, f.data());
    for (double q : f)
        if (!std::isfinite(q))
            throw std::domain_error("Nonfinite vector field");
    return f;
}
PhaseState field(const ODESystem &system, const PhaseState &x) { return fieldAt(system, x, 0, false); }
namespace {
PhaseState derivative(const ODESystem &s, const PhaseState &x, double time, bool variational, bool periodic) {
    auto f = fieldAt(s, x, time, periodic);
    if (!variational)
        return f;
    const auto n = s.descriptor().stateDimension;
    PhaseState a;
    if (periodic) {
        // Differentiate the actual forced RHS, including state-dependent input gains.
        // Some native/extension Jacobians describe only the autonomous drift.
        a.assign(n * n, 0);
        PhaseState plus(x.begin(), x.begin() + n), minus = plus;
        for (std::size_t j = 0; j < n; ++j) {
            const double h = std::cbrt(std::numeric_limits<double>::epsilon()) * std::max(1., std::abs(x[j]));
            plus[j] += h; minus[j] -= h;
            const auto fp = fieldAt(s, plus, time, true), fm = fieldAt(s, minus, time, true);
            for (std::size_t i = 0; i < n; ++i) a[i * n + j] = (fp[i] - fm[i]) / (2 * h);
            plus[j] = minus[j] = x[j];
        }
    }
    else a = s.jacobian(0, x.data());
    if (a.size() != n * n)
        throw std::runtime_error("System Jacobian shape mismatch");
    f.resize(n + n * n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            for (std::size_t k = 0; k < n; ++k)
                f[n + i * n + j] += a[i * n + k] * x[n + k * n + j];
    return f;
}
PhaseState rk4(const ODESystem &s, const PhaseState &x, double time, double h, bool variational, bool periodic) {
    auto a = derivative(s, x, time, variational, periodic), z = x;
    for (std::size_t i = 0; i < x.size(); ++i)
        z[i] = x[i] + h * a[i] / 2;
    auto b = derivative(s, z, time + h / 2, variational, periodic);
    for (std::size_t i = 0; i < x.size(); ++i)
        z[i] = x[i] + h * b[i] / 2;
    auto c = derivative(s, z, time + h / 2, variational, periodic);
    for (std::size_t i = 0; i < x.size(); ++i)
        z[i] = x[i] + h * c[i];
    auto d = derivative(s, z, time + h, variational, periodic);
    for (std::size_t i = 0; i < x.size(); ++i)
        z[i] = x[i] + h * (a[i] + 2 * b[i] + 2 * c[i] + d[i]) / 6;
    return z;
}
} // namespace
Flow flow(const ODESystem &system, const PhaseState &initial, double duration,
          const PhaseSpaceOptions &o, bool variational, bool stopAtDomain,
          const ExecutionContext *execution, double startTime) {
    Flow result;
    const bool periodic = o.dynamics == "periodic-forcing";
    const auto n = initial.size();
    auto x = initial;
    if (variational) {
        x.resize(n + n * n);
        for (std::size_t i = 0; i < n; ++i)
            x[n + i * n + i] = 1;
    }
    result.curve.time.push_back(startTime);
    result.curve.values = initial;
    double t = 0, h = std::min(o.maximumStep, o.outputStep), target = o.outputStep;
    const double end = std::abs(duration), direction = duration < 0 ? -1. : 1.;
    for (std::size_t iteration = 0; t < end && iteration < o.maximumIntegrationSteps; ++iteration) {
        if (execution)
            execution->checkpoint();
        h = std::min({h, end - t, target - t});
        if (h < 1e-14 * (1 + t)) {
            result.reason = "Minimum integration step reached";
            break;
        }
        PhaseState coarse, fine;
        double error = 0;
        try {
            const double now = startTime + direction * t;
            coarse = rk4(system, x, now, direction * h, variational, periodic);
            fine = rk4(system, rk4(system, x, now, direction * h / 2, variational, periodic),
                       now + direction * h / 2, direction * h / 2, variational, periodic);
            for (std::size_t i = 0; i < x.size(); ++i) {
                const double scale =
                    o.absoluteTolerance +
                    o.relativeTolerance * std::max(std::abs(x[i]), std::abs(fine[i]));
                error = std::max(error, std::abs(fine[i] - coarse[i]) / (15 * scale));
            }
        } catch (const std::exception &e) {
            h *= .25;
            result.reason = e.what();
            continue;
        }
        if (!std::isfinite(error)) {
            h *= .25;
            continue;
        }
        if (error > 1) {
            h *= std::max(.1, .9 * std::pow(error, -.2));
            continue;
        }
        for (std::size_t i = 0; i < x.size(); ++i)
            fine[i] += (fine[i] - coarse[i]) / 15;
        if (stopAtDomain && !inside(fine, o)) {
            result.reason = "Reached declared phase-space boundary";
            break;
        }
        x = std::move(fine);
        t += h;
        if (t >= target - 1e-12 * (1 + target) || t >= end - 1e-12 * (1 + end)) {
            result.curve.time.push_back(startTime + direction * t);
            result.curve.values.insert(result.curve.values.end(), x.begin(), x.begin() + n);
            target = std::min(end, target + o.outputStep);
        }
        h = std::min(o.maximumStep,
                     h * std::clamp(error > 0 ? .9 * std::pow(error, -.2) : 2., .5, 2.));
        if (end - t < 1e-13 * (1 + end)) {
            t = end;
            break;
        }
    }
    result.complete = t >= end;
    result.state.assign(x.begin(), x.begin() + n);
    if (variational)
        result.monodromy.assign(x.begin() + n, x.end());
    if (result.complete)
        result.reason.clear();
    else if (result.reason.empty())
        result.reason = "Integration step limit reached";
    result.curve.status = result.complete ? "complete" : "truncated";
    result.curve.reason = result.reason;
    return result;
}
} // namespace dynsyslab::core::phase_detail
