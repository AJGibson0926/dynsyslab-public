#include "PhaseSurfaceInternal.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace dynsyslab::core::phase_detail::surface {
namespace {
struct Shot {
    PhaseState state, parameterDerivative;
    double foot{}, residual{}, speed{};
};
Shot onPlane(const ODESystem &system, const Front &front, const PhaseState &reference,
             const PhaseState &normal, double time, double foot, double sign,
             const PhaseSpaceOptions &o, Budget &budget, const ExecutionContext *execution) {
    const auto n = reference.size();
    for (std::size_t iteration = 0; iteration < o.maximumNewtonIterations; ++iteration) {
        if (execution)
            execution->checkpoint();
        if (++budget.evaluations > budget.maximum)
            throw std::runtime_error("Surface shooting evaluation budget reached");
        const auto initial = interpolate(front.points, foot);
        auto run = flow(system, initial, sign * time, o, true, true, execution);
        if (!run.complete)
            throw std::runtime_error(run.reason + " at surface flight time=" +
                                     std::to_string(time) + ", foot=" + std::to_string(foot));
        if (run.state.size() != n || run.monodromy.size() != n * n)
            throw std::runtime_error("Surface variational flow shape mismatch");
        const double residual = dot(subtract(run.state, reference), normal, o);
        const auto left = interpolate(front.points, foot - 1e-4),
                   right = interpolate(front.points, foot + 1e-4);
        PhaseState mapped(n);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j)
                mapped[i] += run.monodromy[i * n + j] * (right[j] - left[j]) / (2e-4);
        const double derivative = dot(mapped, normal, o);
        if (!std::isfinite(derivative) || std::abs(derivative) < 1e-13)
            throw std::runtime_error("Surface boundary-value chart lost rank");
        if (std::abs(residual) <= o.residualTolerance) {
            const double speed = -sign * dot(field(system, run.state), normal, o) / derivative;
            return {std::move(run.state), std::move(mapped), foot, std::abs(residual), speed};
        }
        double correction = std::clamp(-residual / derivative, -2., 2.);
        bool accepted = false;
        for (int trial = 0; trial < 10; ++trial, correction *= .5) {
            if (++budget.evaluations > budget.maximum)
                throw std::runtime_error("Surface shooting evaluation budget reached");
            const auto candidate = flow(system, interpolate(front.points, foot + correction),
                                        sign * time, o, false, true, execution);
            if (candidate.complete && std::abs(dot(subtract(candidate.state, reference), normal,
                                                   o)) < std::abs(residual)) {
                foot += correction;
                accepted = true;
                break;
            }
        }
        if (!accepted)
            throw std::runtime_error("Surface boundary-value correction did not descend");
    }
    throw std::runtime_error("Surface boundary-value correction did not converge");
}
} // namespace
Point advance(const ODESystem &system, const Front &front, double parent, double distance,
              double sign, const PhaseSpaceOptions &o, Budget &budget,
              const ExecutionContext *execution) {
    const auto reference = interpolate(front.points, parent),
               normal = tangent(front.points, parent, o);
    const auto outward = unit(subtract(reference, interpolate(front.previous, parent)), o);
    auto current = onPlane(system, front, reference, normal, 0, parent, sign, o, budget, execution);
    double time = 0, previousDistance = 0;
    // Continue the trivial orbit and bracket its first outward distance crossing.
    for (int step = 0; step < 512; ++step) {
        auto velocity = field(system, current.state);
        for (std::size_t j = 0; j < velocity.size(); ++j)
            velocity[j] = sign * velocity[j] + current.speed * current.parameterDerivative[j];
        const double outwardSpeed = dot(velocity, outward, o);
        if (outwardSpeed <= 1e-12 * std::max(1., length(velocity, o)))
            throw std::runtime_error("Surface front no longer has an outward regular continuation");
        double dt = .75 * distance / std::max(outwardSpeed, 1e-14);
        dt = std::min(dt, 2. / std::max(1e-6, std::abs(current.speed)));
        dt = std::min(dt, o.horizon - time);
        if (!(dt > 1e-12 * (1 + time)))
            throw std::runtime_error("Surface orbit horizon reached before the next front");
        auto next = onPlane(system, front, reference, normal, time + dt,
                            current.foot + dt * current.speed, sign, o, budget, execution);
        const double nextDistance = length(subtract(next.state, reference), o);
        if (nextDistance >= distance) {
            double low = time, high = time + dt, bestTime = high;
            auto best = next;
            for (int iteration = 0; iteration < 50; ++iteration) {
                const double error = length(subtract(best.state, reference), o) - distance;
                if (std::abs(error) <= o.residualTolerance) {
                    if (dot(subtract(best.state, reference), outward, o) <= 0)
                        throw std::runtime_error("Surface shooting selected an inward crossing");
                    return {best.state, parent, std::max(best.residual, std::abs(error)), bestTime,
                            best.foot};
                }
                if (error >= 0)
                    high = bestTime;
                else
                    low = bestTime;
                auto distanceVelocity = field(system, best.state);
                for (std::size_t j = 0; j < distanceVelocity.size(); ++j)
                    distanceVelocity[j] = sign * distanceVelocity[j] + best.speed * best.parameterDerivative[j];
                const auto displacement = subtract(best.state, reference);
                const double slope =
                    dot(displacement, distanceVelocity, o) / std::max(length(displacement, o), 1e-14);
                // Safeguarded Newton on the distance condition avoids dozens of full orbit solves.
                double mid = bestTime - error / slope;
                if (!std::isfinite(mid) || slope <= 0 || mid <= low || mid >= high)
                    mid = (low + high) / 2;
                bestTime = mid;
                best = onPlane(system, front, reference, normal, mid,
                               current.foot + (mid - time) * current.speed, sign, o, budget,
                               execution);
            }
            throw std::runtime_error("Surface front distance condition did not converge");
        }
        if (nextDistance + o.residualTolerance < previousDistance)
            throw std::runtime_error("Surface continuation folded before the next geodesic front");
        previousDistance = nextDistance;
        time += dt;
        current = std::move(next);
    }
    throw std::runtime_error("Surface continuation step budget reached");
}
} // namespace dynsyslab::core::phase_detail::surface
