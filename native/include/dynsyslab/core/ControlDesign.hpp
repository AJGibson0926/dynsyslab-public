#pragma once
// Classical state-feedback and observer design on an LTI realization:
// pole placement, discrete LQR/LQE (DARE), continuous LQR/LQE (CARE, reusing
// the suite's Hamiltonian solver) and LQG output-feedback compensators.
// KLQR/KLQE (LQR.hpp, Estimation.hpp) remain the Koopman-specific runtime
// controllers; these routines operate on any StateSpaceModel.
#include <dynsyslab/core/StateSpace.hpp>
#include <dynsyslab/numerics/Control.hpp>

#include <complex>
#include <cstddef>
#include <vector>

namespace dynsyslab::core {

// K (inputs x states) with eig(A - B K) = poles. Complex poles must appear in
// conjugate pairs. Single input: Ackermann (repeated poles allowed). Multiple
// inputs: robust eigenvector assignment (Kautsky-Nichols-Van Dooren method 0);
// a pole may repeat at most `inputs` times.
[[nodiscard]] std::vector<double> placePoles(const StateSpaceModel& model,
    const std::vector<std::complex<double>>& poles);
// Observer gain L (states x outputs) with eig(A - L C) = poles (dual problem).
[[nodiscard]] std::vector<double> placeObserverPoles(const StateSpaceModel& model,
    const std::vector<std::complex<double>>& poles);

// The one LQR result: Riccati solution P (states x states), gain K (inputs x states, u = -K x),
// closed-loop poles and the Riccati residual.
// Continuous: CARE. Discrete: DARE by the structure-preserving doubling
// algorithm, K = (R + B^T P B)^-1 B^T P A.
[[nodiscard]] numerics::ContinuousLqrResult designLQR(const StateSpaceModel& model,
    const std::vector<double>& statePenalty, const std::vector<double>& inputPenalty);

struct EstimatorDesign {
    std::vector<double> covariance; // states x states (discrete: a-priori steady state)
    std::vector<double> gain;       // states x outputs: predictor gain (x_hat' = A x_hat + B u + L (y - C x_hat))
    std::vector<double> filterGain; // discrete only: measurement-update gain P C^T (C P C^T + V)^-1
    std::vector<std::complex<double>> observerPoles; // eig(A - L C)
    double relativeResidual{};
};
// W: process covariance (states x states, disturbance enters every state);
// V: measurement covariance (outputs x outputs).
[[nodiscard]] EstimatorDesign designLQE(const StateSpaceModel& model,
    const std::vector<double>& processCovariance, const std::vector<double>& measurementCovariance);

struct LqgDesign {
    numerics::ContinuousLqrResult regulator;
    EstimatorDesign estimator;
    StateSpaceModel compensator; // input y, output u: x_hat' = (A - BK - LC + LDK) x_hat + L y, u = -K x_hat
    StateSpaceModel closedLoop;  // state [x; x_hat], input r added to u, output y
    std::vector<std::complex<double>> closedLoopPoles; // = regulator poles U observer poles (separation)
};
[[nodiscard]] LqgDesign designLQG(const StateSpaceModel& model,
    const std::vector<double>& statePenalty, const std::vector<double>& inputPenalty,
    const std::vector<double>& processCovariance, const std::vector<double>& measurementCovariance);

// Closed loop of plant and static state feedback u = -K x + r.
[[nodiscard]] StateSpaceModel stateFeedbackClosedLoop(const StateSpaceModel& model,
    const std::vector<double>& gain);

}  // namespace dynsyslab::core
