// Prehistory generation and clock comparison (TrajectoryAudit.hpp). Ported from the side-cars
// (miscellaneous_computes prehistory.backward_history; trajectory_resolution compare_clocks).
#include <dynsyslab/core/TrajectoryAudit.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

// RK4 of the unforced field from t over T (either sign) in equal steps <= h (Hairer, Norsett and
// Wanner, Solving ODEs I, 1993, section II.1).
void rk4(const ODESystem& system, std::vector<double>& x, double t, double T, double h) {
    const auto n = x.size();
    const auto steps = std::max<std::size_t>(1, std::size_t(std::ceil(std::abs(T) / h)));
    const double dt = T / double(steps);
    std::vector<double> k1(n), k2(n), k3(n), k4(n), y(n);
    for (std::size_t s = 0; s < steps; ++s) {
        const double at = t + double(s) * dt;
        system.rhsWithInput(at, x.data(), 0, k1.data());
        for (std::size_t i = 0; i < n; ++i) y[i] = x[i] + dt / 2 * k1[i];
        system.rhsWithInput(at + dt / 2, y.data(), 0, k2.data());
        for (std::size_t i = 0; i < n; ++i) y[i] = x[i] + dt / 2 * k2[i];
        system.rhsWithInput(at + dt / 2, y.data(), 0, k3.data());
        for (std::size_t i = 0; i < n; ++i) y[i] = x[i] + dt * k3[i];
        system.rhsWithInput(at + dt, y.data(), 0, k4.data());
        for (std::size_t i = 0; i < n; ++i) x[i] += dt / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
        for (const auto v : x) if (!std::isfinite(v)) throw std::runtime_error("prehistory: the backward flow diverged");
    }
}

}  // namespace

Prehistory backwardHistory(const ODESystem& system, const std::vector<double>& x0, double t0, const std::vector<double>& times,
                           double maximumStep) {
    const auto n = system.descriptor().stateDimension;
    if (x0.size() != n || times.empty() || !(maximumStep > 0)) throw std::invalid_argument("prehistory: state, times and step");
    for (std::size_t k = 0; k < times.size(); ++k)
        if (!std::isfinite(times[k]) || (k && !(times[k] > times[k - 1])) || !(times[k] < t0))
            throw std::invalid_argument("prehistory times must increase and lie strictly before the first sample");
    Prehistory out;
    out.states.resize(times.size() * n);
    auto x = x0;
    double at = t0;
    for (std::size_t k = times.size(); k-- > 0;) {
        rk4(system, x, at, times[k] - at, maximumStep);
        at = times[k];
        std::copy(x.begin(), x.end(), out.states.begin() + std::ptrdiff_t(k * n));
    }
    rk4(system, x, at, t0 - at, maximumStep);  // forward check
    for (std::size_t i = 0; i < n; ++i) out.roundtripError = std::max(out.roundtripError, std::abs(x[i] - x0[i]));
    return out;
}

ClockComparison compareClocks(const std::vector<double>& rt, const std::vector<double>& ref, const std::vector<double>& st,
                              const std::vector<double>& sampled, std::size_t d, double begin, double end) {
    const auto increasing = [](const std::vector<double>& t) {
        if (t.size() < 2) return false;
        for (std::size_t k = 1; k < t.size(); ++k) if (!(t[k] > t[k - 1])) return false;
        return true;
    };
    if (!d || !increasing(rt) || !increasing(st) || ref.size() != rt.size() * d || sampled.size() != st.size() * d)
        throw std::invalid_argument("clock comparison needs increasing clocks and matching sample-by-coordinate arrays");
    ClockComparison out;
    out.first = std::max({rt.front(), st.front(), begin});
    out.last = std::min({rt.back(), st.back(), end});
    const auto lo = std::size_t(std::lower_bound(rt.begin(), rt.end(), out.first) - rt.begin());
    const auto hi = std::size_t(std::upper_bound(rt.begin(), rt.end(), out.last) - rt.begin());
    if (lo >= hi) throw std::invalid_argument("clock comparison: no reference samples in the common window");
    out.referenceCount = hi - lo;
    out.sampleCount = std::size_t(std::upper_bound(st.begin(), st.end(), out.last) - std::lower_bound(st.begin(), st.end(), out.first));
    out.maximumAbsolute.assign(d, 0); out.rms.assign(d, 0); out.worstTime.assign(d, rt[lo]);
    for (auto k = lo; k < hi; ++k) {
        const auto j = std::min(std::size_t(std::upper_bound(st.begin(), st.end(), rt[k]) - st.begin()), st.size() - 1);
        const auto a = j ? j - 1 : 0;
        const double f = st[j] == st[a] ? 0 : (rt[k] - st[a]) / (st[j] - st[a]);
        for (std::size_t c = 0; c < d; ++c) {
            const double e = sampled[a * d + c] + f * (sampled[j * d + c] - sampled[a * d + c]) - ref[k * d + c];
            out.rms[c] += e * e;
            if (std::abs(e) > out.maximumAbsolute[c]) { out.maximumAbsolute[c] = std::abs(e); out.worstTime[c] = rt[k]; }
        }
    }
    for (auto& r : out.rms) r = std::sqrt(r / double(out.referenceCount));
    return out;
}

}  // namespace dynsyslab::core
