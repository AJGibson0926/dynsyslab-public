#pragma once
#include <dynsyslab/core/ShootingIdentification.hpp>
namespace dynsyslab::core {
struct SequentialShootingOptions {
    // Indices in the uniformly sampled observation stream. Each update uses
    // exactly the prefix [0,index]. At least two samples must be available.
    std::vector<std::size_t> updateIndices;
    std::size_t predictionSubsteps{1};
    bool retainOriginalStarts{false};
};
struct SequentialShootingUpdate {
    std::size_t measurementIndex{};
    double time{};
    ShootingIdentificationResult identification;
};
struct SequentialShootingResult {
    std::vector<double> time,states,informationTime;
    std::size_t stateDimension{};
    std::vector<SequentialShootingUpdate> updates;
};
// Expanding-window deterministic estimation. Refit initial parameters with
// each available prefix, then predict only until the next scheduled update.
// Emitted samples are never revised using later observations. This is not a
// Kalman filter, finite-window MHE, uncertainty estimator, or real-time claim.
// Uses single-shooting weighted least squares (Nocedal & Wright, 2006, ch.10).
SequentialShootingResult sequentialShooting(ODESystem&,
    const std::vector<double>& times,const std::vector<double>& measurements,
    const std::vector<double>& variances,const std::vector<ShootingStart>&,
    const InitialStateMap&,const InputAwareReadout&,const KnownScalarInput&,
    const SequentialShootingOptions&,IntegrationOptions = {},
    const numerics::NonlinearLeastSquaresOptions& = {});
}
