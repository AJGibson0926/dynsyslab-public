#pragma once
// Balanced-truncation sidecar for Koopman control (KLQR).
//
// Balances (A, B, C_w) where C_w is the performance output (identity unless
// measurements are named), designs LQR on the reduced model, and lifts the
// gain back as K_full = K_r P so it can drive the full lifted state. The
// full-model closed-loop poles eig(A - B K_r P) are the acceptance check.
#include <dynsyslab/core/BalancedTruncation.hpp>
#include <dynsyslab/core/ControlDesign.hpp>
#include <dynsyslab/core/KoopmanStateSpace.hpp>

#include <complex>
#include <vector>

namespace dynsyslab::core {

struct ReducedRegulatorDesign {
    BalancedReduction reduction;
    numerics::ContinuousLqrResult regulator;      // on the reduced model
    std::vector<double> liftedGain; // inputs x full states: K_r P
    std::vector<std::complex<double>> liftedClosedLoopPoles; // eig(A - B K_r P)
    bool liftedClosedLoopStable{};
};

// State penalty in FULL coordinates (n x n), projected as L^T Q L; input
// penalty inputs x inputs.
[[nodiscard]] ReducedRegulatorDesign designReducedLQR(const StateSpaceModel& model,
    const BalancedReductionOptions& reduction,
    const std::vector<double>& statePenalty, const std::vector<double>& inputPenalty);

// KROM convenience: options.measurements empty -> performance output C_w = I.
[[nodiscard]] ReducedRegulatorDesign designReducedKoopmanLQR(const KoopmanModel& model,
    const KoopmanStateSpaceOptions& realization, const BalancedReductionOptions& reduction,
    const std::vector<double>& statePenalty, const std::vector<double>& inputPenalty);

}  // namespace dynsyslab::core
