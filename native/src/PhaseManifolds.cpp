#include "PhaseSpaceInternal.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
namespace dynsyslab::core::phase_detail {
namespace {
double difference(const PhaseState &a, const PhaseState &b, const PhaseSpaceOptions &o) {
    auto d = a;
    for (std::size_t j = 0; j < d.size(); ++j)
        d[j] -= b[j];
    return norm(d, o);
}
double geometricError(const PhaseCurve &a, const PhaseCurve &b, const PhaseSpaceOptions &o) {
    const auto n = o.lower.size();
    if (a.time.size() < 2 || b.time.size() < 2)
        return std::numeric_limits<double>::infinity();
    double maximum = 0;
    // Compare the shared branch extent; a truncated refinement does not certify its tail.
    for (std::size_t i = 0; i < a.time.size(); ++i) {
        double best = std::numeric_limits<double>::infinity();
        for (std::size_t k = 1; k < b.time.size(); ++k) {
            double dot = 0, length = 0;
            for (std::size_t j = 0; j < n; ++j) {
                double scale = o.upper[j] - o.lower[j],
                       v = (b.values[k * n + j] - b.values[(k - 1) * n + j]) / scale;
                dot += (a.values[i * n + j] - b.values[(k - 1) * n + j]) / scale * v;
                length += v * v;
            }
            double u = length > 0 ? std::clamp(dot / length, 0., 1.) : 0;
            PhaseState p(n), q(n);
            for (std::size_t j = 0; j < n; ++j) {
                p[j] = a.values[i * n + j];
                q[j] = b.values[(k - 1) * n + j] +
                       u * (b.values[k * n + j] - b.values[(k - 1) * n + j]);
            }
            best = std::min(best, difference(p, q, o));
        }
        if (std::isfinite(best))
            maximum = std::max(maximum, best);
    }
    return maximum;
}
} // namespace
void traceManifolds(const ODESystem &system, PhaseSpaceResult &result,
                    const ExecutionContext *execution) {
    const auto &o = result.options;
    const auto n = system.descriptor().stateDimension;
    for (std::size_t index = 0; index < result.fixedPoints.size(); ++index) {
        const auto &point = result.fixedPoints[index];
        if (point.stability != "saddle")
            continue;
        for (bool stable : {true, false}) {
            const auto count = stable ? point.stableDimension : point.unstableDimension;
            if (count == 2 && o.surface.enabled) {
                traceSurface(system, result, index, stable, execution);
                continue;
            }
            if (count != 1) {
                result.diagnostics.push_back("Equilibrium " + std::to_string(index) +
                                             (stable ? " stable" : " unstable") +
                                             " manifold has dimension " + std::to_string(count) +
                                             "; surface continuation is disabled or requires a different manifold dimension");
                continue;
            }
            std::size_t column = 0;
            for (; column < n; ++column)
                if (stable ? point.eigenvalues[column].real() < -o.eigenvalueTolerance
                           : point.eigenvalues[column].real() > o.eigenvalueTolerance)
                    break;
            if (column == n || std::abs(point.eigenvalues[column].imag()) > o.eigenvalueTolerance)
                continue;
            PhaseState tangent(n);
            for (std::size_t j = 0; j < n; ++j)
                tangent[j] = point.eigenvectors[j * n + column].real();
            const double size = norm(tangent, o);
            if (size == 0)
                continue;
            for (auto &q : tangent)
                q /= size;
            const double direction = stable ? -1. : 1.;
            for (double sign : {-1., 1.}) {
                auto initial = point.state, refined = initial;
                for (std::size_t j = 0; j < n; ++j) {
                    initial[j] += sign * o.seedOffset * tangent[j];
                    refined[j] += sign * o.seedOffset * .5 * tangent[j];
                }
                if (!inside(initial, o) || !inside(refined, o))
                    continue;
                auto trace =
                    flow(system, initial, direction * o.horizon, o, false, true, execution);
                const double warmup = std::log(2.) / std::abs(point.eigenvalues[column].real());
                if (warmup > o.horizon) {
                    trace.curve.reason += "; refinement timescale exceeds horizon";
                    trace.curve.status = "unresolved";
                } else {
                    auto warm =
                        flow(system, refined, direction * warmup, o, false, true, execution);
                    if (warm.complete) {
                        auto reference = flow(system, warm.state, direction * o.horizon, o, false,
                                              true, execution);
                        trace.curve.refinementError =
                            geometricError(trace.curve, reference.curve, o);
                        if (trace.curve.refinementError > 10 * o.seedOffset) {
                            trace.curve.status = "unresolved";
                            trace.curve.reason += "; seed-halving geometry did not converge";
                        }
                    } else {
                        trace.curve.status = "unresolved";
                        trace.curve.reason +=
                            "; seed refinement could not reach comparison section";
                    }
                }
                trace.curve.kind = stable ? "stable_manifold" : "unstable_manifold";
                trace.curve.baseObject = index;
                trace.curve.seedOffset = o.seedOffset;
                // Local central-difference tangency check on the sampled curve.
                for (std::size_t k = 1; k + 1 < trace.curve.time.size(); ++k) {
                    PhaseState x(n), velocity(n);
                    double dt = trace.curve.time[k + 1] - trace.curve.time[k - 1];
                    for (std::size_t j = 0; j < n; ++j) {
                        x[j] = trace.curve.values[k * n + j];
                        velocity[j] = (trace.curve.values[(k + 1) * n + j] -
                                       trace.curve.values[(k - 1) * n + j]) /
                                      dt;
                    }
                    try {
                        auto f = field(system, x);
                        trace.curve.tangentDefect =
                            std::max(trace.curve.tangentDefect,
                                     difference(velocity, f, o) / (1 + norm(f, o)));
                    } catch (const std::exception &) {
                        trace.curve.status = "unresolved";
                        trace.curve.reason += "; tangent check encountered invalid dynamics";
                        break;
                    }
                }
                result.manifolds.push_back(std::move(trace.curve));
            }
        }
    }
}
} // namespace dynsyslab::core::phase_detail
