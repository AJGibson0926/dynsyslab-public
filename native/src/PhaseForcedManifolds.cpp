#include "PhaseSpaceInternal.hpp"
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
namespace dynsyslab::core::phase_detail {
namespace {
using Line = std::vector<PhaseState>;
PhaseState difference(PhaseState a, const PhaseState &b) {
    for (std::size_t j = 0; j < a.size(); ++j)
        a[j] -= b[j];
    return a;
}
PhaseState midpoint(const PhaseState &a, const PhaseState &b) {
    auto x = a;
    for (std::size_t j = 0; j < x.size(); ++j)
        x[j] = (a[j] + b[j]) / 2;
    return x;
}
struct Branch {
    Line points;
    std::string status{"complete"}, reason;
    double defect{};
};
Branch trace(const ODESystem &system, const PhasePeriodicOrbit &orbit, const PhaseState &direction,
             double offset, double duration, std::size_t iterations, const PhaseSpaceOptions &o,
             const ExecutionContext *execution) {
    Branch branch;
    std::size_t evaluations = 0;
    auto map = [&](const PhaseState &x) {
        if (++evaluations > 8 * o.maximumManifoldPoints)
            throw std::runtime_error("Return-map evaluation budget reached");
        auto run = flow(system, x, duration, o, false, true, execution, o.forcingSectionTime);
        if (!run.complete)
            throw std::runtime_error(run.reason);
        return run.state;
    };
    try {
        auto first = orbit.sectionState;
        for (std::size_t j = 0; j < first.size(); ++j)
            first[j] += offset * direction[j];
        if (!inside(first, o))
            throw std::runtime_error("Local manifold seed leaves the declared domain");
        Line segment{first, map(first)};
        branch.points = {orbit.sectionState, first, segment.back()};
        for (std::size_t step = 0; step < iterations; ++step) {
            if (execution)
                execution->checkpoint();
            Line next;
            const auto firstImage = map(segment.front());
            next.push_back(firstImage);
            std::function<void(const PhaseState &, const PhaseState &, const PhaseState &,
                               const PhaseState &, int)>
                refine;
            refine = [&](const PhaseState &a, const PhaseState &b, const PhaseState &fa,
                         const PhaseState &fb, int depth) {
                if (next.size() + branch.points.size() >= o.maximumManifoldPoints)
                    throw std::runtime_error("Manifold point budget reached");
                const auto middle = midpoint(a, b), fm = map(middle);
                const double defect = norm(difference(fm, midpoint(fa, fb)), o);
                const double spacing = norm(difference(fb, fa), o);
                if (defect > o.manifoldTolerance || spacing > std::sqrt(o.manifoldTolerance)) {
                    if (depth >= 18)
                        throw std::runtime_error(
                            "Manifold interpolation refinement did not converge");
                    refine(a, middle, fa, fm, depth + 1);
                    refine(middle, b, fm, fb, depth + 1);
                } else {
                    branch.defect = std::max(branch.defect, defect);
                    next.push_back(fb);
                }
            };
            for (std::size_t j = 1; j < segment.size(); ++j) {
                const auto fa = next.back(), fb = map(segment[j]);
                refine(segment[j - 1], segment[j], fa, fb, 0);
            }
            branch.points.insert(branch.points.end(), next.begin() + 1, next.end());
            segment = std::move(next);
        }
    } catch (const ExecutionCancelled &) {
        throw;
    } catch (const std::exception &e) {
        branch.status = "truncated";
        branch.reason = e.what();
    }
    return branch;
}
double geometricError(const Line &coarse, const Line &fine, const PhaseSpaceOptions &o) {
    if (coarse.size() < 3 || fine.size() < 3)
        return std::numeric_limits<double>::infinity();
    double error = 0;
    for (const auto &x : coarse) {
        double best = std::numeric_limits<double>::infinity();
        for (std::size_t j = 1; j < fine.size(); ++j) {
            auto delta = difference(fine[j], fine[j - 1]), gap = difference(x, fine[j - 1]);
            double dot = 0, length = 0;
            for (std::size_t k = 0; k < x.size(); ++k) {
                const double width = o.upper[k] - o.lower[k];
                dot += gap[k] * delta[k] / (width * width);
                length += delta[k] * delta[k] / (width * width);
            }
            const double u = length > 0 ? std::clamp(dot / length, 0., 1.) : 0;
            for (std::size_t k = 0; k < x.size(); ++k)
                gap[k] -= u * delta[k];
            best = std::min(best, norm(gap, o));
        }
        error = std::max(error, best);
    }
    return error;
}
} // namespace
void traceForcedManifolds(const ODESystem &system, PhaseSpaceResult &result,
                          const ExecutionContext *execution) {
    const auto &o = result.options;
    const auto n = system.descriptor().stateDimension;
    for (std::size_t index = 0; index < result.periodicOrbits.size(); ++index) {
        const auto &orbit = result.periodicOrbits[index];
        if (orbit.stability != "saddle_periodic_response")
            continue;
        for (bool stable : {true, false}) {
            const auto dimension = stable ? orbit.stableDimension : orbit.unstableDimension;
            if (dimension != 1) {
                result.diagnostics.push_back(
                    "Periodic-response manifold has dimension greater than one on the forcing "
                    "section; surface continuation is required");
                continue;
            }
            for (std::size_t k = 0; k < n; ++k) {
                const auto lambda = orbit.floquetMultipliers[k];
                if (std::abs(lambda.imag()) > o.eigenvalueTolerance ||
                    (stable ? std::abs(lambda) >= 1 - 10 * o.eigenvalueTolerance
                            : std::abs(lambda) <= 1 + 10 * o.eigenvalueTolerance))
                    continue;
                PhaseState direction(n);
                for (std::size_t j = 0; j < n; ++j)
                    direction[j] = orbit.floquetVectors[j * n + k].real();
                const double scale = norm(direction, o);
                if (scale == 0)
                    continue;
                for (auto &v : direction)
                    v /= scale;
                // Negative multipliers exchange branches; the squared map preserves their
                // orientation.
                const double turns = lambda.real() < 0 ? 2. : 1.;
                const double duration = (stable ? -1. : 1.) * turns * orbit.period;
                const double growth =
                    std::pow(stable ? 1 / std::abs(lambda) : std::abs(lambda), turns);
                const auto warmup = std::min<std::size_t>(
                    64, std::size_t(std::ceil(std::log(2.) / std::log(growth))) + 2);
                for (double sign : {-1., 1.}) {
                    auto coarse = trace(system, orbit, direction, sign * o.seedOffset, duration,
                                        o.manifoldIterations, o, execution);
                    auto tighter = o;
                    tighter.manifoldTolerance *= .25;
                    tighter.maximumStep *= .5;
                    tighter.relativeTolerance *= .125;
                    tighter.absoluteTolerance *= .125;
                    auto fine = trace(system, orbit, direction, sign * o.seedOffset / 2, duration,
                                      o.manifoldIterations + warmup, tighter, execution);
                    PhaseCurve curve;
                    curve.kind = stable ? "stable_section_manifold" : "unstable_section_manifold";
                    curve.baseKind = "periodic_response";
                    curve.baseObject = index;
                    curve.seedOffset = o.seedOffset;
                    curve.tangentDefect = coarse.defect;
                    curve.refinementError = geometricError(coarse.points, fine.points, o);
                    curve.status = curve.refinementError <= 10 * o.manifoldTolerance ? coarse.status
                                                                                     : "unresolved";
                    curve.reason = curve.status == "unresolved"
                                       ? "Forcing-section curve failed seed-halving and "
                                         "interpolation refinement"
                                       : coarse.reason;
                    // These samples parameterize a curve in one forcing section, not a time
                    // trajectory.
                    for (const auto &x : coarse.points)
                        curve.values.insert(curve.values.end(), x.begin(), x.end());
                    result.manifolds.push_back(std::move(curve));
                }
            }
        }
    }
}
} // namespace dynsyslab::core::phase_detail
