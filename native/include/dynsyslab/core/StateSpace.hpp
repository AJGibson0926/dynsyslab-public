#pragma once
// Classical linear time-invariant state-space toolkit.
//
//   continuous:  x' = A x + B u,        y = C x + D u
//   discrete:    x[k+1] = A x[k] + B u[k], y[k] = C x[k] + D u[k]
//
// Matrices are row-major std::vector<double>. Every Koopman model DynSysLab
// produces (K-ROM split realization, H-DMD, DMDc/EDMDc, linear SINDy part) is
// an ordinary LTI realization once a measurement C and an input map B are
// chosen; KoopmanStateSpace.hpp performs that conversion.
#include <complex>
#include <cstddef>
#include <string>
#include <vector>

namespace dynsyslab::core {

enum class TimeDomain { continuous, discrete };

struct StateSpaceModel {
    std::vector<double> A; // states x states
    std::vector<double> B; // states x inputs (may be empty when inputs == 0)
    std::vector<double> C; // outputs x states (may be empty when outputs == 0)
    std::vector<double> D; // outputs x inputs
    std::size_t states{}, inputs{}, outputs{};
    TimeDomain domain{TimeDomain::continuous};
    double sampleTime{}; // > 0 for discrete models
};

// Validates dimensions and finiteness. An empty D is replaced by zeros.
[[nodiscard]] StateSpaceModel makeStateSpace(
    std::vector<double> A, std::vector<double> B, std::vector<double> C,
    std::vector<double> D, std::size_t states, std::size_t inputs,
    std::size_t outputs, TimeDomain domain = TimeDomain::continuous,
    double sampleTime = 0.0);
void validateStateSpace(const StateSpaceModel& model);
[[nodiscard]] const char* timeDomainName(TimeDomain domain) noexcept;

[[nodiscard]] std::vector<std::complex<double>> systemPoles(const StateSpaceModel& model);
// Continuous: max Re(lambda) < -margin. Discrete: max |mu| < 1 - margin.
[[nodiscard]] bool isAsymptoticallyStable(const StateSpaceModel& model, double margin = 0.0);

// --- Lyapunov / Stein equations (Bartels-Stewart on the complex Schur form) ---
// Continuous: A X + X A^T + Q = 0.  Requires lambda_i + conj(lambda_j) != 0.
[[nodiscard]] std::vector<double> solveContinuousLyapunov(
    const std::vector<double>& A, const std::vector<double>& Q, std::size_t states);
// Discrete (Stein): A X A^T - X + Q = 0.  Requires mu_i conj(mu_j) != 1.
[[nodiscard]] std::vector<double> solveDiscreteLyapunov(
    const std::vector<double>& A, const std::vector<double>& Q, std::size_t states);

// --- Gramians ---------------------------------------------------------------
// Infinite-horizon gramians require an asymptotically stable model.
//   Wc: A Wc + Wc A^T + B B^T = 0  (discrete: A Wc A^T - Wc + B B^T = 0)
//   Wo: A^T Wo + Wo A + C^T C = 0  (discrete: A^T Wo A - Wo + C^T C = 0)
[[nodiscard]] std::vector<double> controllabilityGramian(const StateSpaceModel& model);
[[nodiscard]] std::vector<double> observabilityGramian(const StateSpaceModel& model);
// Finite horizon: continuous integral over [0, horizon] (Van Loan block
// exponential, valid for unstable models); discrete sum of `horizon` terms
// k = 0 .. horizon-1 (horizon must be a positive integer).
[[nodiscard]] std::vector<double> finiteControllabilityGramian(
    const StateSpaceModel& model, double horizon);
[[nodiscard]] std::vector<double> finiteObservabilityGramian(
    const StateSpaceModel& model, double horizon);

// Empirical (simulation/data) gramians in the sense of Lall, Marsden and
// Glavaski: impulse responses of every input (controllability) and output
// differences of +/- initial-state perturbations (observability), integrated
// with the trapezoid rule on `samples` points of [0, horizon]. For an LTI
// model they converge to the finite-horizon gramians at O(dt^2); for a
// nonlinear plant callers supply recorded responses to snapshotGramian.
struct EmpiricalGramianOptions {
    double horizon{};          // continuous time span, or discrete step count
    std::size_t samples{1001}; // continuous quadrature points (ignored for discrete)
    std::vector<double> scales{1.0}; // perturbation magnitudes, averaged
};
[[nodiscard]] std::vector<double> empiricalControllabilityGramian(
    const StateSpaceModel& model, const EmpiricalGramianOptions& options);
[[nodiscard]] std::vector<double> empiricalObservabilityGramian(
    const StateSpaceModel& model, const EmpiricalGramianOptions& options);
// Data-driven gramian sum_k w_k x_k x_k^T from recorded snapshots
// (dimension x samples, row-major). Empty weights means unit weights.
[[nodiscard]] std::vector<double> snapshotGramian(
    const std::vector<double>& snapshots, std::size_t dimension,
    std::size_t samples, const std::vector<double>& weights = {});

// --- Discretization and simulation -------------------------------------------
enum class DiscretizationMethod { zeroOrderHold, tustin, forwardEuler };
[[nodiscard]] StateSpaceModel discretize(const StateSpaceModel& continuous,
    double sampleTime, DiscretizationMethod method = DiscretizationMethod::zeroOrderHold);
// Inverse Tustin (bilinear) map; exact inverse of discretize(..., tustin).
[[nodiscard]] StateSpaceModel continuousFromTustin(const StateSpaceModel& discrete);

struct LinearResponse {
    std::vector<double> time;
    std::vector<double> outputs; // samples x outputs
    std::vector<double> states;  // samples x states
    std::size_t samples{};
};
// Inputs are samples x inputs, zero-order held between uniform samples.
[[nodiscard]] LinearResponse simulateLinear(const StateSpaceModel& model,
    const std::vector<double>& time, const std::vector<double>& inputs,
    const std::vector<double>& initialState = {});
[[nodiscard]] LinearResponse impulseResponse(const StateSpaceModel& model,
    std::size_t input, double horizon, std::size_t samples);
[[nodiscard]] LinearResponse stepResponse(const StateSpaceModel& model,
    std::size_t input, double horizon, std::size_t samples);

// --- Stable / unstable additive decomposition ---------------------------------
// Real spectral projection with the matrix sign function: G = G_s + G_u, with
// x = T [x_s; x_u]. Modes with Re(lambda) >= -margin (|mu| >= 1 - margin) are
// placed in the unstable part, so marginal Koopman modes are retained intact.
struct StableUnstableSplit {
    StateSpaceModel stable, unstable; // D belongs to the stable part
    std::vector<double> transform, inverse; // states x states
    std::size_t stableStates{}, unstableStates{};
    double margin{};
};
[[nodiscard]] StableUnstableSplit splitStableUnstable(
    const StateSpaceModel& model, double margin = 1.0e-8);

// Parallel connection G1 + sign * G2 (sign = -1 gives the error system G1 - G2).
// A model with zero states contributes only its D term.
[[nodiscard]] StateSpaceModel parallelConnection(const StateSpaceModel& first,
    const StateSpaceModel& second, double sign = 1.0);

// Similarity transform x = T z: (T^-1 A T, T^-1 B, C T, D).
[[nodiscard]] StateSpaceModel similarityTransform(const StateSpaceModel& model,
    const std::vector<double>& transform, const std::vector<double>& inverse);

}  // namespace dynsyslab::core
