#pragma once
// Classical robustness of SISO feedback loops: loops broken at the plant input or
// output for state feedback (LQR) and output feedback (LQG compensators), classical
// and disk margins, sensitivity peaks, and loop transfer recovery (Doyle-Stein).
//
// Doyle (1978): with a scalar input the LQR loop K (sI - A)^-1 B satisfies
// |1 + L(jw)| >= 1 (gain margin [1/2, inf), phase margin >= 60 deg); an LQG loop,
// in which a Kalman filter replaces the state, has no guaranteed margin.
#include <dynsyslab/core/ControlDesign.hpp>
#include <dynsyslab/core/SystemNorms.hpp>

#include <complex>
#include <vector>

namespace dynsyslab::core {

struct LoopRobustness {
    // Classical margins from the Nyquist curve's crossings of the negative real axis:
    // the loop tolerates gains in (gainMarginLower, gainMarginUpper) about nominal 1.
    double gainMarginUpper{};     // smallest crossing factor above 1 (infinite if none)
    double gainMarginLower{};     // largest crossing factor below 1 (0 if none)
    double phaseMargin{};         // degrees in [0, 180]: angular distance from -1 at the nearest |L| = 1 crossing (infinite if none)
    double phaseCrossover{}, gainCrossover{};  // rad per model time unit
    // Disk (return difference) measures: alpha = min |1 + L(jw)| = 1 / ||S||_inf.
    double minimumReturnDifference{};
    double returnDifferenceFrequency{};
    double diskGainLower{}, diskGainUpper{};  // [1/(1+alpha), 1/(1-alpha)] (upper infinite if alpha >= 1)
    double diskPhaseMargin{};                 // degrees, 2 asin(min(alpha, 2) / 2)
    double sensitivityPeak{};                 // ||1 / (1 + L)||_inf
    double complementarySensitivityPeak{};    // ||L / (1 + L)||_inf
    bool closedLoopStable{};                  // nominal loop
    std::vector<double> frequencies;          // evaluation grid
    std::vector<std::complex<double>> response;  // L(jw) on that grid
};

// L(s) = K (sI - A)^-1 B: state feedback u = -K x, broken at the (single) input.
[[nodiscard]] StateSpaceModel stateFeedbackInputLoop(const StateSpaceModel& plant, const std::vector<double>& gain);

// L(s) = C (sI - A)^-1 L_e: the Kalman filter loop at the plant's single output (dual of
// the LQR loop; with scalar V it obeys the same |1 + L| >= 1 guarantee).
[[nodiscard]] StateSpaceModel estimatorOutputLoop(const StateSpaceModel& plant, const std::vector<double>& estimatorGain);

// L(s) = -C_K(s) G(s) for a compensator u = C_K(s) y (e.g. LqgDesign::compensator),
// broken at the plant's single input; states [x; x_c].
[[nodiscard]] StateSpaceModel outputFeedbackInputLoop(const StateSpaceModel& plant, const StateSpaceModel& compensator);

// L(s) = -G(s) C_K(s), broken at the plant's single output; states [x; x_c].
[[nodiscard]] StateSpaceModel outputFeedbackOutputLoop(const StateSpaceModel& plant, const StateSpaceModel& compensator);

// Classical, disk and sensitivity measures of a SISO negative-feedback loop.
[[nodiscard]] LoopRobustness loopRobustness(const StateSpaceModel& loop, std::size_t gridPoints = 4000);

// Loop transfer recovery at the plant input (Doyle-Stein): W + q2 B B^T as the
// estimator's process covariance. As q2 grows, the LQG input loop approaches the
// LQR loop for minimum-phase plants (check transmissionZeros).
[[nodiscard]] std::vector<double> loopTransferRecoveryCovariance(const StateSpaceModel& plant,
    const std::vector<double>& processCovariance, double q2);

}  // namespace dynsyslab::core
