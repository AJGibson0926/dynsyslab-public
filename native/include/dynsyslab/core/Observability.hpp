#pragma once
#include <cstddef>
#include <dynsyslab/core/Pbh.hpp>
#include <vector>
#include <string>
namespace dynsyslab::core {
class KoopmanModel;
struct ReconstructionResult;
enum class ObservabilityMethod {
    classicalMatrix,
    discreteGramian,
    discreteInfiniteGramian,
    continuousFiniteGramian,
    continuousInfiniteGramian
};

struct ObservabilityOptions {
    ObservabilityMethod method{ObservabilityMethod::classicalMatrix};
    std::size_t blockCount{};
    double horizon{1.0};
    double sampleTime{};
    double relativeTolerance{1.0e-10};
    double maximumAbsoluteEntry{1.0e150};
};

struct ObservabilityResult {
    std::vector<double> matrix;              // block rows x state
    std::vector<double> gramian;             // state x state
    std::vector<double> singularValues; // sqrt eigenvalues of Gramian
    std::vector<double> gramianEigenvalues;
    std::vector<double> leftDirections; // matrix rows x min(rows,states), thin SVD
    std::vector<double> rightDirections;     // state x state, columns ordered by observability
    ObservabilityMethod method{ObservabilityMethod::classicalMatrix};
    std::size_t stateDimension{};
    std::size_t outputDimension{};
    std::size_t blockCount{};
    std::size_t rank{};
    double rankTolerance{};
    double horizon{};
    double conditionNumber{}; // infinite when rank deficient
    double resolvedConditionNumber{}; // largest / smallest retained strength
    double effectiveRank{};
};

[[nodiscard]] ObservabilityResult analyzeObservability(
    const std::vector<double>& generator, const std::vector<double>& output,
    std::size_t stateDimension, std::size_t outputDimension,
    std::size_t blockCount = 0, double relativeTolerance = 1.0e-10);

[[nodiscard]] ObservabilityResult analyzeObservability(
    const std::vector<double>& generator, const std::vector<double>& output,
    std::size_t stateDimension, std::size_t outputDimension,
    const ObservabilityOptions& options);

[[nodiscard]] const char* observabilityMethodName(ObservabilityMethod method) noexcept;


// KROM continuous generator by default. For discrete methods sampleTime must
// be positive; it specifies the transition's dt and is applied exactly once.
[[nodiscard]] ObservabilityResult analyzeObservability(const KoopmanModel& model,
    const std::string& reconstruction, const ObservabilityOptions& options = {});

// Fractional squared participation of coordinate groups in each observability direction:
// entry [g, j] = sum over i in group g of V[i, j]^2, with V the orthonormal right directions
// (state x state, row-major), so every column sums to one over a partition of the coordinates.
// Groups must be disjoint and inside the state. Returns groups x state, row-major.
[[nodiscard]] std::vector<double> observabilityModeParticipation(
    const ObservabilityResult& result, const std::vector<std::vector<std::size_t>>& groups);
// The KROM's split coordinates grouped by nonredundant Koopman mode: a real mode is one
// coordinate, a complex mode its real and imaginary coordinates (KoopmanModel::splitValues order).
[[nodiscard]] std::vector<std::vector<std::size_t>> splitModeGroups(const KoopmanModel& model);

// Convert G = C Psi into a real output map over the KROM's split lifted state.
// Defined in dynsyslab_linear_analysis; shared by observability, the
// linear-systems toolkit and KLQE.
[[nodiscard]] std::vector<double> reconstructionSplitMatrix(
    const KoopmanModel& model, const ReconstructionResult& reconstruction);
} // namespace dynsyslab::core
