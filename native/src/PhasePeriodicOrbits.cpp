#include "PhaseSpaceInternal.hpp"
#include <algorithm>
#include <cmath>
#include <dynsyslab/numerics/Eigensystem.hpp>
#include <dynsyslab/numerics/LeastSquares.hpp>
#include <stdexcept>
namespace dynsyslab::core::phase_detail {
namespace {
PhaseState difference(PhaseState x, const PhaseState &y) {
    for (std::size_t i = 0; i < x.size(); ++i)
        x[i] -= y[i];
    return x;
}
std::vector<PeriodicOrbitSeed> candidates(const ODESystem &system, const PhaseSpaceOptions &o,
                                          const ExecutionContext *execution) {
    auto result = o.periodicSeeds;
    if (!result.empty())
        return result;
    auto initial = seeds(o, system.descriptor().stateDimension);
    const auto n = system.descriptor().stateDimension;
    const auto section = sectionGeometry(o, n);
    // Oriented section crossings supply guesses only; shooting validates every orbit.
    for (std::size_t k = 0; k < std::min<std::size_t>(12, initial.size()); ++k) {
        auto warm = flow(system, initial[k], o.transientTime, o, false, true, execution);
        if (!warm.complete)
            continue;
        auto run = flow(system, warm.state, o.horizon, o, false, true, execution);
        double last = -1;
        PhaseState lastState;
        for (std::size_t i = 1; i < run.curve.time.size(); ++i) {
            PhaseState before(run.curve.values.begin() + (i - 1) * n,
                              run.curve.values.begin() + i * n);
            PhaseState after(run.curve.values.begin() + i * n,
                             run.curve.values.begin() + (i + 1) * n);
            const double a = sectionValue(before, section), b = sectionValue(after, section);
            const bool crossing = o.sectionDirection > 0 ? a < 0 && b >= 0
                                : o.sectionDirection < 0 ? a > 0 && b <= 0
                                                         : a * b <= 0 && a != b;
            if (!crossing)
                continue;
            const double u = -a / (b - a), time = run.curve.time[i - 1] +
                                                  u * (run.curve.time[i] - run.curve.time[i - 1]);
            PhaseState state(n);
            for (std::size_t j = 0; j < n; ++j)
                state[j] = run.curve.values[(i - 1) * n + j] +
                           u * (run.curve.values[i * n + j] - run.curve.values[(i - 1) * n + j]);
            if (last >= 0 && time - last > 2 * o.outputStep) {
                result.push_back({lastState, time - last});
                break;
            }
            last = time;
            lastState = std::move(state);
        }
    }
    return result;
}
PhasePeriodicOrbit shoot(const ODESystem &system, const PeriodicOrbitSeed &seed,
                         const PhaseSpaceOptions &o, const ExecutionContext *execution) {
    auto x = seed.state;
    double period = seed.period;
    const auto n = x.size();
    const auto section = sectionGeometry(o, n);
    double sectionScale = 0;
    for (double value : section.scaledNormal)
        sectionScale += value * value;
    sectionScale = std::sqrt(sectionScale);
    auto phase = [&](const PhaseState &q) {
        return sectionValue(q, section) / sectionScale;
    };
    Flow run;
    bool converged = false;
    for (std::size_t iteration = 0; iteration < o.maximumNewtonIterations; ++iteration) {
        if (execution)
            execution->checkpoint();
        run = flow(system, x, period, o, true, false, execution);
        if (!run.complete)
            throw std::runtime_error(run.reason);
        const auto residual = difference(run.state, x);
        const double cost = std::max(norm(residual, o), std::abs(phase(x)));
        if (cost < o.residualTolerance) {
            converged = true;
            break;
        }
        auto f = field(system, run.state);
        std::vector<double> matrix((n + 1) * (n + 1)), rhs(n + 1);
        for (std::size_t i = 0; i < n; ++i) {
            const double width = o.upper[i] - o.lower[i];
            rhs[i] = -residual[i] / width;
            for (std::size_t j = 0; j < n; ++j)
                matrix[i * (n + 1) + j] = (run.monodromy[i * n + j] - (i == j ? 1. : 0.)) *
                                          (o.upper[j] - o.lower[j]) / width;
            matrix[i * (n + 1) + n] = f[i] / width;
            matrix[n * (n + 1) + i] = section.scaledNormal[i] / sectionScale;
        }
        rhs[n] = -phase(x);
        const auto correction =
            dynsyslab::numerics::solveLeastSquares(matrix, n + 1, n + 1, rhs, 1);
        bool accepted = false;
        double alpha = 1;
        for (int trial = 0; trial < 12; ++trial, alpha *= .5) {
            auto next = x;
            for (std::size_t j = 0; j < n; ++j)
                next[j] += alpha * correction[j] * (o.upper[j] - o.lower[j]);
            const double nextPeriod = period + alpha * correction[n];
            if (nextPeriod <= 2 * o.outputStep || nextPeriod > 2 * o.horizon || !inside(next, o))
                continue;
            auto candidate = flow(system, next, nextPeriod, o, false, false, execution);
            if (candidate.complete && std::max(norm(difference(candidate.state, next), o),
                                               std::abs(phase(next))) < cost) {
                x = std::move(next);
                period = nextPeriod;
                accepted = true;
                break;
            }
        }
        if (!accepted)
            break;
    }
    if (!converged)
        throw std::runtime_error("Phase-conditioned shooting did not converge");
    if (norm(field(system, x), o) < 100 * o.residualTolerance)
        throw std::runtime_error("Shooting converged to an equilibrium");
    for (std::size_t i = 0; i < run.curve.time.size(); ++i) {
        PhaseState q(run.curve.values.begin() + i * n, run.curve.values.begin() + (i + 1) * n);
        if (!inside(q, o))
            throw std::runtime_error("Periodic orbit leaves the declared domain");
    }
    auto refinedOptions = o;
    refinedOptions.maximumStep *= .5;
    refinedOptions.relativeTolerance *= .125;
    refinedOptions.absoluteTolerance *= .125;
    auto refined = flow(system, x, period, refinedOptions, true, false, execution);
    if (!refined.complete)
        throw std::runtime_error("Refined autonomous return failed: " + refined.reason);
    const double refinementError = norm(difference(refined.state, run.state), o);
    double variationalError = 0;
    for (std::size_t j = 0; j < run.monodromy.size(); ++j)
        variationalError = std::max(
            variationalError, std::abs(run.monodromy[j] - refined.monodromy[j]) /
                                  (1 + std::abs(refined.monodromy[j])));
    if (refinementError > 20 * o.residualTolerance || variationalError > 1e-5)
        throw std::runtime_error(
            "Autonomous orbit failed tighter-integration return or variational validation");
    const auto returnMap = poincareFromTransition(system, x, refined, o, section);
    if (returnMap.status != "complete")
        throw std::runtime_error(returnMap.reason);
    if (o.sectionDirection && returnMap.direction != o.sectionDirection)
        throw std::runtime_error("Periodic orbit converged on the opposite oriented section branch");
    PhasePeriodicOrbit orbit;
    orbit.period = period;
    orbit.sectionState = x;
    orbit.curve = std::move(refined.curve);
    orbit.curve.kind = "periodic_orbit";
    orbit.curve.refinementError = refinementError;
    orbit.closureResidual = norm(difference(refined.state, x), o);
    orbit.phaseResidual = std::abs(phase(x));
    orbit.variationalRefinementError = variationalError;
    orbit.floquetMultipliers =
        dynsyslab::numerics::eigensystem(refined.monodromy, n, false).eigenvalues;
    std::size_t neutral = 0;
    for (std::size_t j = 1; j < n; ++j)
        if (std::abs(orbit.floquetMultipliers[j] - 1.) <
            std::abs(orbit.floquetMultipliers[neutral] - 1.))
            neutral = j;
    orbit.neutralMultiplierError = std::abs(orbit.floquetMultipliers[neutral] - 1.);
    const auto initialField = field(system, x), finalField = field(system, refined.state);
    PhaseState tangentDefect(n, 0);
    for (std::size_t row = 0; row < n; ++row) {
        tangentDefect[row] = -finalField[row];
        for (std::size_t column = 0; column < n; ++column)
            tangentDefect[row] += refined.monodromy[row * n + column] * initialField[column];
    }
    orbit.tangentIdentityError = norm(tangentDefect, o);
    if (orbit.tangentIdentityError > 1e-5)
        throw std::runtime_error("Autonomous tangent identity M f(x)=f(P(x)) failed validation");
    orbit.returnTime = returnMap.returnTime;
    orbit.returnResidual = returnMap.fixedPointResidual;
    orbit.transversality = returnMap.transversality;
    orbit.transversalityAngle = returnMap.transversalityAngle;
    orbit.returnDerivative = returnMap.returnDerivative;
    orbit.reducedReturnDerivative = returnMap.reducedDerivative;
    orbit.sectionTangentBasis = returnMap.tangentBasis;
    orbit.sectionMultipliers = returnMap.multipliers;
    std::size_t stable = 0, unstable = 0, center = 0;
    for (const auto multiplier : orbit.sectionMultipliers) {
        const double radius = std::abs(multiplier);
        if (radius < 1 - 10 * o.eigenvalueTolerance)
            ++stable;
        else if (radius > 1 + 10 * o.eigenvalueTolerance)
            ++unstable;
        else
            ++center;
    }
    orbit.stableDimension = stable; orbit.unstableDimension = unstable; orbit.centerDimension = center;
    orbit.curve.baseKind = "periodic_orbit";
    orbit.isolated = center == 0;
    orbit.stability = center              ? "nonhyperbolic_periodic_orbit"
                      : stable == orbit.sectionMultipliers.size()   ? "attracting_limit_cycle"
                      : unstable == orbit.sectionMultipliers.size() ? "repelling_limit_cycle"
                                          : "saddle_limit_cycle";
    return orbit;
}
} // namespace
void findPeriodicOrbits(const ODESystem &system, PhaseSpaceResult &result,
                        const ExecutionContext *execution) {
    const auto &o = result.options;
    const auto n = system.descriptor().stateDimension;
    const auto guesses = candidates(system, o, execution);
    for (const auto &seed : guesses) {
        ++result.periodicSearches;
        try {
            auto orbit = shoot(system, seed, o, execution);
            bool duplicate = false;
            for (const auto &previous : result.periodicOrbits) {
                if (std::abs(previous.period - orbit.period) > 1e-4 * orbit.period)
                    continue;
                PhaseState first(orbit.curve.values.begin(), orbit.curve.values.begin() + n);
                for (std::size_t i = 0; i < previous.curve.time.size(); ++i) {
                    PhaseState q(previous.curve.values.begin() + i * n,
                                 previous.curve.values.begin() + (i + 1) * n);
                    if (norm(difference(first, q), o) < .01) {
                        duplicate = true;
                        break;
                    }
                }
                if (duplicate)
                    break;
            }
            if (!duplicate)
                result.periodicOrbits.push_back(std::move(orbit));
        } catch (const ExecutionCancelled &) {
            throw;
        } catch (const std::exception &e) {
            ++result.periodicFailures;
            result.diagnostics.push_back(std::string("Periodic seed: ") + e.what());
        }
    }
    if (guesses.empty())
        result.diagnostics.push_back("No oriented section returns found; this does not establish "
                                     "absence of periodic orbits");
}
} // namespace dynsyslab::core::phase_detail
