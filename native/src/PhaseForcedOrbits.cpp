#include "PhaseSpaceInternal.hpp"
#include <algorithm>
#include <cmath>
#include <dynsyslab/numerics/Eigensystem.hpp>
#include <dynsyslab/numerics/LeastSquares.hpp>
#include <stdexcept>
namespace dynsyslab::core::phase_detail {
namespace {
PhaseState difference(PhaseState a, const PhaseState &b) {
    for (std::size_t j = 0; j < a.size(); ++j)
        a[j] -= b[j];
    return a;
}
PhasePeriodicOrbit shoot(const ODESystem &system, PhaseState x, const PhaseSpaceOptions &o,
                         const ExecutionContext *execution) {
    const auto n = x.size();
    double period = o.forcingPeriod * o.forcingPeriods;
    Flow run;
    bool converged = false;
    for (std::size_t iteration = 0; iteration < o.maximumNewtonIterations; ++iteration) {
        run = flow(system, x, period, o, true, true, execution, o.forcingSectionTime);
        if (!run.complete)
            throw std::runtime_error(run.reason);
        auto residual = difference(run.state, x);
        const double cost = norm(residual, o);
        if (cost < o.residualTolerance) {
            converged = true;
            break;
        }
        std::vector<double> matrix(n * n), rhs(n);
        for (std::size_t i = 0; i < n; ++i) {
            const double width = o.upper[i] - o.lower[i];
            rhs[i] = -residual[i] / width;
            for (std::size_t j = 0; j < n; ++j)
                matrix[i * n + j] = (run.monodromy[i * n + j] - (i == j ? 1. : 0.)) *
                                    (o.upper[j] - o.lower[j]) / width;
        }
        const auto correction = dynsyslab::numerics::solveLeastSquares(matrix, n, n, rhs, 1);
        bool accepted = false;
        double alpha = 1;
        for (int trial = 0; trial < 12; ++trial, alpha *= .5) {
            auto next = x;
            for (std::size_t j = 0; j < n; ++j)
                next[j] += alpha * correction[j] * (o.upper[j] - o.lower[j]);
            if (!inside(next, o))
                continue;
            auto candidate =
                flow(system, next, period, o, false, true, execution, o.forcingSectionTime);
            if (candidate.complete && norm(difference(candidate.state, next), o) < cost) {
                x = std::move(next);
                accepted = true;
                break;
            }
        }
        if (!accepted)
            break;
    }
    if (!converged)
        throw std::runtime_error("Fixed-forcing-phase shooting did not converge");
    // Identify lower-period responses rather than reporting every requested multiple as minimal.
    for (std::size_t divisor = 1; divisor < o.forcingPeriods; ++divisor) {
        if (o.forcingPeriods % divisor)
            continue;
        auto candidate = flow(system, x, o.forcingPeriod * divisor, o, true, true, execution,
                              o.forcingSectionTime);
        if (candidate.complete && norm(difference(candidate.state, x), o) < o.residualTolerance) {
            period = o.forcingPeriod * divisor;
            run = std::move(candidate);
            break;
        }
    }
    auto refinedOptions = o;
    refinedOptions.maximumStep *= .5;
    refinedOptions.relativeTolerance *= .125;
    refinedOptions.absoluteTolerance *= .125;
    auto refined =
        flow(system, x, period, refinedOptions, true, true, execution, o.forcingSectionTime);
    if (!refined.complete)
        throw std::runtime_error("Refined forcing return failed: " + refined.reason);
    const double refinementError = norm(difference(refined.state, run.state), o);
    const double closure = norm(difference(refined.state, x), o);
    if (refinementError > 10 * o.residualTolerance || closure > 20 * o.residualTolerance)
        throw std::runtime_error("Periodic response failed tighter-integration return validation");
    double variationalError = 0;
    for (std::size_t j = 0; j < run.monodromy.size(); ++j)
        variationalError =
            std::max(variationalError, std::abs(run.monodromy[j] - refined.monodromy[j]) /
                                           (1 + std::abs(refined.monodromy[j])));
    if (!std::isfinite(variationalError) || variationalError > 1e-5)
        throw std::runtime_error("Floquet transition matrix failed tighter-integration validation");
    PhasePeriodicOrbit orbit;
    orbit.period = period;
    orbit.returnTime = period;
    orbit.sectionState = x;
    orbit.curve = std::move(refined.curve);
    orbit.curve.kind = "forced_periodic_response";
    orbit.curve.baseKind = "periodic_response";
    orbit.curve.refinementError = refinementError;
    orbit.closureResidual = closure;
    orbit.returnResidual = closure;
    orbit.variationalRefinementError = variationalError;
    const auto spectral = dynsyslab::numerics::eigensystem(refined.monodromy, n, true);
    orbit.floquetMultipliers = spectral.eigenvalues;
    orbit.floquetVectors = spectral.rightEigenvectors;
    orbit.sectionMultipliers = spectral.eigenvalues;
    orbit.returnDerivative = refined.monodromy;
    orbit.reducedReturnDerivative = refined.monodromy;
    orbit.sectionTangentBasis.assign(n * n, 0);
    for (std::size_t j = 0; j < n; ++j)
        orbit.sectionTangentBasis[j * n + j] = 1;
    // All multipliers are transverse on the forcing section; no autonomous neutral multiplier is
    // removed.
    for (auto value : orbit.floquetMultipliers) {
        if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
            throw std::runtime_error("Nonfinite forcing-section Floquet multiplier");
        const double radius = std::abs(value);
        if (radius < 1 - 10 * o.eigenvalueTolerance)
            ++orbit.stableDimension;
        else if (radius > 1 + 10 * o.eigenvalueTolerance)
            ++orbit.unstableDimension;
        else
            ++orbit.centerDimension;
    }
    orbit.isolated = orbit.centerDimension == 0;
    orbit.stability = orbit.centerDimension          ? "nonhyperbolic_periodic_response"
                      : orbit.stableDimension == n   ? "attracting_periodic_response"
                      : orbit.unstableDimension == n ? "repelling_periodic_response"
                                                     : "saddle_periodic_response";
    return orbit;
}
} // namespace
void findForcedPeriodicOrbits(const ODESystem &system, PhaseSpaceResult &result,
                              const ExecutionContext *execution) {
    const auto &o = result.options;
    auto guesses = std::vector<PhaseState>{};
    for (const auto &seed : o.periodicSeeds)
        guesses.push_back(seed.state);
    if (guesses.empty()) {
        guesses = seeds(o, system.descriptor().stateDimension);
        if (guesses.size() > 32) {
            std::vector<PhaseState> sampled;
            for (std::size_t j = 0; j < 32; ++j)
                sampled.push_back(guesses[j * (guesses.size() - 1) / 31]);
            guesses = std::move(sampled);
        }
    }
    for (const auto &seed : guesses) {
        ++result.periodicSearches;
        try {
            auto orbit = shoot(system, seed, o, execution);
            bool duplicate = false;
            for (const auto &previous : result.periodicOrbits)
                if (norm(difference(previous.sectionState, orbit.sectionState), o) <
                        100 * o.residualTolerance &&
                    std::abs(previous.period - orbit.period) < o.residualTolerance * orbit.period)
                    duplicate = true;
            if (!duplicate)
                result.periodicOrbits.push_back(std::move(orbit));
        } catch (const ExecutionCancelled &) {
            throw;
        } catch (const std::exception &e) {
            ++result.periodicFailures;
            result.diagnostics.push_back(std::string("Forcing-section seed: ") + e.what());
        }
    }
    if (result.periodicOrbits.empty())
        result.diagnostics.push_back("No periodic response was verified in this finite search; "
                                     "this does not establish absence");
}
} // namespace dynsyslab::core::phase_detail
