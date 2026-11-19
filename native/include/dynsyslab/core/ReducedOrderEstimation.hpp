#pragma once
// Balanced-truncation sidecar for Koopman estimation (KLQE).
//
// Large H-DMD/echo KROMs carry many weakly observable modes. The sidecar
// balances the measurement model (A, G, C) where G is the process-noise input
// (identity unless a named control map is used), designs a steady-state LQE
// on the reduced model, and reports what the reduced gain does when lifted
// back into the full model. The established full-order KLQE is unchanged.
#include <dynsyslab/core/BalancedTruncation.hpp>
#include <dynsyslab/core/ControlDesign.hpp>
#include <dynsyslab/core/Estimation.hpp>
#include <dynsyslab/core/KoopmanStateSpace.hpp>

#include <complex>
#include <vector>

namespace dynsyslab::core {

struct ReducedEstimatorDesign {
    BalancedReduction reduction;
    EstimatorDesign estimator;    // on the reduced model
    std::vector<double> liftedGain; // full states x outputs: lift * L_r
    std::vector<std::complex<double>> liftedObserverPoles; // eig(A - lift L_r C) of the full model
    bool liftedObserverStable{};
};

// Process covariance is given in FULL coordinates (n x n) and projected as
// P W P^T; measurement covariance is outputs x outputs.
[[nodiscard]] ReducedEstimatorDesign designReducedLQE(const StateSpaceModel& model,
    const BalancedReductionOptions& reduction,
    const std::vector<double>& processCovariance, const std::vector<double>& measurementCovariance);

// KROM convenience: options.controlMap empty -> balancing input G = I.
[[nodiscard]] ReducedEstimatorDesign designReducedKoopmanLQE(const KoopmanModel& model,
    const KoopmanStateSpaceOptions& realization, const BalancedReductionOptions& reduction,
    const std::vector<double>& processCovariance, const std::vector<double>& measurementCovariance);

// KLQE stage option `balanced_truncation`: balance the design's (A, I, C_y),
// design the LQE on the reduced model and write the lifted result back into
// `design.solution` (gain L L_r, observer poles of the full KROM, covariance
// L P_r L^T). The KLQE runtime then runs unchanged at full order. Throws when the
// lifted observer is not asymptotically stable.
ReducedEstimatorDesign reduceKoopmanLQEDesign(KoopmanLQEDesign& design,
    const BalancedReductionOptions& reduction,
    const std::vector<double>& processCovariance, const std::vector<double>& measurementCovariance);

}  // namespace dynsyslab::core
