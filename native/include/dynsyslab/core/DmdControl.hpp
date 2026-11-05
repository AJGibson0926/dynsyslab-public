#pragma once

#include <complex>
#include <cstddef>
#include <vector>
#include <dynsyslab/core/SINDy.hpp>
#include <dynsyslab/core/PolynomialEigenfunctions.hpp>

namespace dynsyslab::core {

struct DmdControlOptions {
    std::size_t augmentedRank{}; // Rank of [X; U]; zero selects numerical rank.
    std::size_t rank{};          // Output POD rank; zero selects numerical rank.
    double singularValueTolerance{}; // Absolute threshold; zero is scale-aware.
};

// All matrices are row-major. Identification uses q_next = A q + B u.
// q denotes supplied measurements, NOT necessarily physical state.
struct DmdControlResult {
    std::size_t dimension{}, inputDimension{}, pairs{}, rank{}, augmentedRank{};
    std::size_t augmentedNumericalRank{};
    double deltaTime{}, oneStepRelativeError{}, projectionRelativeError{};
    double retainedConditionNumber{}, eigenvectorConditionNumber{};
    bool fullRowRank{}; // Necessary, not sufficient, for physical identifiability.
    bool hasEigenfunctions{}; // False if the reduced eigenbasis is ill-conditioned.
    std::vector<double> basis;          // dimension x rank: z = P^T q
    std::vector<double> operatorMatrix; // rank x rank: reduced discrete A
    std::vector<double> controlMatrix;  // rank x inputDimension: reduced discrete B
    std::vector<double> augmentedSingularValues, outputSingularValues;
    std::vector<std::complex<double>> discreteEigenvalues;
    std::vector<std::complex<double>> modes; // dimension x rank: P W
    std::vector<std::complex<double>> exactModes; // dimension x rank: A_hat P W
    // Rows L = W^-1 P^T, NOT right eigenvectors. psi(q) = L q.
    // These are eigenfunctions of the fitted reduced model, not certified
    // Koopman eigenfunctions of the original (possibly non-Markov) measurements.
    std::vector<std::complex<double>> eigenfunctionCoefficients; // rank x dimension
    std::vector<std::complex<double>> modalInputMatrix; // rank x inputDimension
    // Export linear polynomials to the shared eigenfunction representation.
    // Uses principal log(mu)/dt; rejects zero multipliers/defective eigenbases.
    [[nodiscard]] PolynomialEigenfunctionSet eigenfunctionSet() const;
    // Compatibility export requires finite continuous rates (nonzero multipliers).
    [[nodiscard]] KoopmanSINDyResult polynomialEigenfunctions() const;
};

// Explicit sample-major pairs support multiple trajectories without ever
// connecting endpoints. inputs[k] is applied on X[k] -> Y[k], not at Y[k].
// No derivatives, equations of motion, delay construction or centering occur.
[[nodiscard]] DmdControlResult dynamicModeDecompositionWithControl(
    const std::vector<double>& x, const std::vector<double>& y,
    const std::vector<double>& inputs, std::size_t pairs,
    std::size_t dimension, std::size_t inputDimension, double deltaTime,
    const DmdControlOptions& options = {});

// Stateless evaluation in the SAME snapshot coordinates/units used to fit.
// Any echo-history buffer and frozen latent transform belong to the caller.
[[nodiscard]] std::vector<std::complex<double>> evaluateDmdControlEigenfunctions(
    const DmdControlResult& model, const std::vector<double>& snapshots,
    std::size_t samples);
[[nodiscard]] std::vector<double> predictDmdControl(
    const DmdControlResult& model, const std::vector<double>& snapshots,
    const std::vector<double>& inputs, std::size_t samples);

} // namespace dynsyslab::core
