#include "PhaseSpaceInternal.hpp"
#include <algorithm>
#include <cmath>
#include <dynsyslab/numerics/Eigensystem.hpp>
#include <dynsyslab/numerics/LeastSquares.hpp>
namespace dynsyslab::core::phase_detail {
std::vector<PhaseState> seeds(const PhaseSpaceOptions &o, std::size_t n) {
    auto result = o.equilibriumSeeds;
    std::size_t total = 1;
    for (std::size_t j = 0; j < n && total <= o.maximumSeeds; ++j)
        total *= o.gridPerAxis;
    const auto count = std::min(o.maximumSeeds, total);
    std::vector<unsigned> primes;
    for (unsigned candidate = 2; primes.size() < n; ++candidate) {
        bool prime = true;
        for (unsigned p : primes)
            if (candidate % p == 0) {
                prime = false;
                break;
            }
        if (prime)
            primes.push_back(candidate);
    }
    for (std::size_t k = 0; k < count; ++k) {
        PhaseState x(n);
        std::size_t index = k;
        for (std::size_t j = 0; j < n; ++j) {
            double fraction = 0;
            if (total <= o.maximumSeeds) {
                fraction = double(index % o.gridPerAxis) / (o.gridPerAxis - 1);
                index /= o.gridPerAxis;
            } else {
                std::size_t q = k + 1;
                double divisor = 1;
                while (q) {
                    divisor *= primes[j];
                    fraction += (q % primes[j]) / divisor;
                    q /= primes[j];
                }
            }
            x[j] = o.lower[j] + fraction * (o.upper[j] - o.lower[j]);
        }
        result.push_back(std::move(x));
    }
    return result;
}
void findFixedPoints(const ODESystem &system, PhaseSpaceResult &result,
                     const ExecutionContext *execution) {
    const auto &o = result.options;
    const auto n = system.descriptor().stateDimension;
    const auto starts = seeds(o, n);
    for (const auto &start : starts) {
        if (execution)
            execution->report({"phase/fixed-points", "Finding equilibrium candidates",
                               double(result.equilibriumSearches) / starts.size(),
                               result.equilibriumSearches, starts.size()});
        ++result.equilibriumSearches;
        auto x = start;
        bool converged = false;
        try {
            for (std::size_t iteration = 0; iteration < o.maximumNewtonIterations; ++iteration) {
                const auto f = field(system, x);
                const double residual = norm(f, o);
                if (residual <= o.residualTolerance) {
                    converged = true;
                    break;
                }
                auto a = system.jacobian(0, x.data()), b = f;
                for (auto &q : b)
                    q = -q;
                // Solve in domain-scaled coordinates to avoid unit-dependent pivots.
                for (std::size_t i = 0; i < n; ++i) {
                    b[i] /= o.upper[i] - o.lower[i];
                    for (std::size_t j = 0; j < n; ++j)
                        a[i * n + j] *= (o.upper[j] - o.lower[j]) / (o.upper[i] - o.lower[i]);
                }
                const auto delta = dynsyslab::numerics::solveLeastSquares(a, n, n, b, 1);
                double alpha = 1;
                bool accepted = false;
                for (int trial = 0; trial < 18; ++trial, alpha *= .5) {
                    auto next = x;
                    for (std::size_t j = 0; j < n; ++j)
                        next[j] += alpha * delta[j] * (o.upper[j] - o.lower[j]);
                    if (inside(next, o) && norm(field(system, next), o) < residual) {
                        x = std::move(next);
                        accepted = true;
                        break;
                    }
                }
                if (!accepted)
                    break;
            }
            if (!converged)
                converged = norm(field(system, x), o) <= o.residualTolerance;
            if (!converged || !inside(x, o)) {
                ++result.equilibriumFailures;
                continue;
            }
            bool duplicate = false;
            for (const auto &known : result.fixedPoints) {
                auto dx = x;
                for (std::size_t j = 0; j < n; ++j)
                    dx[j] -= known.state[j];
                if (norm(dx, o) < std::max(10 * o.residualTolerance, 1e-6)) {
                    duplicate = true;
                    break;
                }
            }
            if (duplicate)
                continue;
            const auto eigen = dynsyslab::numerics::eigensystem(system.jacobian(0, x.data()), n);
            PhaseFixedPoint point;
            point.state = x;
            point.residual = norm(field(system, x), o);
            point.eigenvalues = eigen.eigenvalues;
            point.eigenvectors = eigen.rightEigenvectors;
            for (auto lambda : point.eigenvalues) {
                if (lambda.real() < -o.eigenvalueTolerance)
                    ++point.stableDimension;
                else if (lambda.real() > o.eigenvalueTolerance)
                    ++point.unstableDimension;
                else
                    ++point.centerDimension;
            }
            point.stability = point.centerDimension          ? "nonhyperbolic"
                              : point.stableDimension == n   ? "sink"
                              : point.unstableDimension == n ? "source"
                                                             : "saddle";
            result.fixedPoints.push_back(std::move(point));
        } catch (const ExecutionCancelled &) {
            throw;
        } catch (const std::exception &) {
            ++result.equilibriumFailures;
        }
    }
}
} // namespace dynsyslab::core::phase_detail
