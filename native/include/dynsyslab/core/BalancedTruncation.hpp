#pragma once
// Hankel singular values, balanced realizations and balanced model reduction,
// plus the data-driven realization tools that produce balanced models from
// recorded responses (ERA, OKID, balanced POD).
#include <dynsyslab/core/StateSpace.hpp>

#include <cstddef>
#include <vector>

namespace dynsyslab::core {

// sigma_i = sqrt(eig(Wc Wo)), descending, from the square-root algorithm.
// Requires an asymptotically stable model.
[[nodiscard]] std::vector<double> hankelSingularValues(const StateSpaceModel& model);

struct BalancedRealization {
    StateSpaceModel balanced;   // order k <= n: numerically minimal balanced part
    std::vector<double> hankelSingularValues; // all n values, descending
    std::vector<double> transform; // n x k, x = T z
    std::vector<double> inverse;   // k x n, z = T^+ x (T^+ T = I)
    std::size_t order{};
};
// Square-root balancing. Directions with sigma <= relativeTolerance*sigma_1
// are uncontrollable or unobservable and are discarded (k < n).
[[nodiscard]] BalancedRealization balancedRealization(
    const StateSpaceModel& model, double relativeTolerance = 1.0e-12);

enum class BalancedReductionMethod { truncation, residualization };
[[nodiscard]] const char* balancedReductionMethodName(BalancedReductionMethod method) noexcept;

struct BalancedReductionOptions {
    // Total reduced order including retained unstable/marginal modes. Zero
    // selects the smallest order meeting errorTolerance (or, if that is also
    // zero, the numerically minimal order).
    std::size_t order{};
    double errorTolerance{};
    BalancedReductionMethod method{BalancedReductionMethod::truncation};
    // Modes with Re(lambda) >= -stabilityMargin (|mu| >= 1 - margin) are kept
    // exactly and excluded from balancing, as is standard for unstable plants.
    double stabilityMargin{1.0e-8};
    bool retainUnstable{true};
    double hankelRelativeTolerance{1.0e-12};
};

struct BalancedReduction {
    StateSpaceModel reduced;
    std::vector<double> hankelSingularValues; // stable part, descending
    std::vector<double> projection; // r x n: z_r = P x
    std::vector<double> lift;       // n x r: x ~ L z_r (P L = I)
    std::size_t order{}, stableOrder{}, unstableOrder{}, fullOrder{};
    double errorBound{}; // ||G - G_r||_inf <= 2 sum of discarded sigma_i
    double lowerBound{}; // ||G - G_r||_inf >= first discarded sigma
    BalancedReductionMethod method{};
    // Stable parts before/after reduction (zero states when absent). Their
    // difference is the whole approximation error, since unstable modes are kept.
    StateSpaceModel stablePart, reducedStablePart;
};
[[nodiscard]] BalancedReduction balancedReduction(
    const StateSpaceModel& model, const BalancedReductionOptions& options = {});

// --- Data-driven realization --------------------------------------------------
// Markov parameters are samples x outputs x inputs, row-major:
// Y_0 = D, Y_k = C A^(k-1) B.
struct EraOptions {
    std::size_t order{};        // 0: rank by relativeTolerance
    std::size_t blockRows{};    // 0: (samples-1)/2
    std::size_t blockColumns{}; // 0: samples-1-blockRows
    double relativeTolerance{1.0e-10};
    double sampleTime{1.0};
};
struct EraResult {
    StateSpaceModel model; // discrete, balanced with respect to the finite Hankel matrix
    std::vector<double> hankelSingularValues;
    std::size_t order{};
};
[[nodiscard]] EraResult eigensystemRealization(const std::vector<double>& markov,
    std::size_t samples, std::size_t outputs, std::size_t inputs, const EraOptions& options = {});

// Observer/Kalman filter identification: Markov parameters from arbitrary
// input/output records (samples x inputs, samples x outputs) through a
// least-squares ARX observer of `observerLag` steps.
[[nodiscard]] std::vector<double> observerKalmanMarkovParameters(
    const std::vector<double>& inputs, const std::vector<double>& outputs,
    std::size_t samples, std::size_t inputCount, std::size_t outputCount,
    std::size_t observerLag, std::size_t markovCount);

// Balanced POD (Rowley 2005): primal snapshots X (n x a, e.g. weighted impulse
// responses) and adjoint snapshots Z (n x b). SVD of Z^T X gives approximate
// Hankel singular values and bi-orthogonal modes with Psi^T Phi = I.
struct BalancedPod {
    std::vector<double> modes;        // n x r (Phi)
    std::vector<double> adjointModes; // n x r (Psi)
    std::vector<double> hankelSingularValues;
    std::size_t states{}, order{};
};
[[nodiscard]] BalancedPod balancedPod(const std::vector<double>& primal,
    const std::vector<double>& adjoint, std::size_t states, std::size_t primalSamples,
    std::size_t adjointSamples, std::size_t order = 0, double relativeTolerance = 1.0e-10);
// Petrov-Galerkin projection (Psi^T A Phi, Psi^T B, C Phi, D) with lift = Phi
// and adjointModes = Psi, both states x order.
[[nodiscard]] StateSpaceModel projectStateSpace(const StateSpaceModel& model,
    const std::vector<double>& lift, const std::vector<double>& adjointModes, std::size_t order);

}  // namespace dynsyslab::core
