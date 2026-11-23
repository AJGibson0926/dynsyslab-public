#pragma once
#include <dynsyslab/core/KnownInput.hpp>
#include <dynsyslab/numerics/NonlinearLeastSquares.hpp>
namespace dynsyslab::core {
using InitialStateMap=std::function<std::vector<double>(const std::vector<double>&)>;
using InputAwareReadout=std::function<std::vector<double>(double,const double*,double)>;
struct ShootingStart { std::vector<double> initial,lower,upper; };
struct ShootingAttempt {
    numerics::NonlinearLeastSquaresResult fit;
    bool valid{};
    std::string failure;
};
struct ShootingIdentificationResult {
    numerics::NonlinearLeastSquaresResult best;
    std::vector<ShootingAttempt> attempts;
    std::size_t bestStart{},trajectoryEvaluations{};
};
// Single shooting over a supplied PAST, uniformly sampled measurement window.
// Unknowns parameterize the initial state; input is known and state-independent.
// Observation may depend directly on input. Independent Gaussian errors give
// variance-weighted least squares. No future data or system initialState() used.
// Local damped Gauss-Newton: Nocedal & Wright (2006), Numerical Optimization, ch.10.
// Failed trajectories invalidate a start, not the entire multi-start run; any
// finite initial candidate is retained if its optimizer fails, with diagnostics.
ShootingIdentificationResult identifyByShooting(ODESystem&,
    const std::vector<double>& times,const std::vector<double>& measurements,
    const std::vector<double>& variances,const std::vector<ShootingStart>&,
    const InitialStateMap&,const InputAwareReadout&,const KnownScalarInput&,
    IntegrationOptions = {},const numerics::NonlinearLeastSquaresOptions& = {});
}
