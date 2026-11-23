#pragma once
#include <functional>
#include <string>
#include <vector>
namespace dynsyslab::numerics {
struct NonlinearLeastSquaresOptions {
    std::size_t maximumIterations{100};
    double finiteDifferenceStep{1e-5},stepTolerance{1e-9},gradientTolerance{1e-8};
};
struct NonlinearLeastSquaresResult {
    std::vector<double> parameters;
    double squaredResidual{};
    std::size_t iterations{},evaluations{};
    bool converged{};
    std::string reason;
};
// Box-bounded damped Gauss-Newton with QR augmented solves, scaled to unit box.
// Finite-difference Jacobian. Local optimization only; no global certificate.
// Nocedal & Wright, Numerical Optimization (2006), chapter 10.
NonlinearLeastSquaresResult boundedNonlinearLeastSquares(
    const std::function<std::vector<double>(const std::vector<double>&)>& residual,
    const std::vector<double>& initial,const std::vector<double>& lower,
    const std::vector<double>& upper,const NonlinearLeastSquaresOptions& = {});
}
