#pragma once
#include <dynsyslab/core/KoopmanModel.hpp>
namespace dynsyslab::core {
enum class EigenfunctionValidationInput { physicalState, coordinates };
struct EigenfunctionErrorStatistics {
    std::size_t samplePairs{};
    double absoluteRms{}, maximumAbsolute{}, relativeL2{};
    double medianPointwiseRelative{}, p95PointwiseRelative{};
};
struct EigenfunctionHorizonError {
    double horizon{};
    EigenfunctionErrorStatistics statistics;
};
struct EigenfunctionModeValidation {
    std::complex<double> eigenvalue; // legacy continuous rate; consult availability
    EigenfunctionErrorStatistics oneStep, anchoredLongTerm, generator;
    std::vector<EigenfunctionHorizonError> horizons;
    double oneStepChangeRelativeL2{};
    bool oneStepChangeDefined{};
    bool continuousRateAvailable{true}, generatorAvailable{true};
};
// Evidence for ONE evaluated dataset, not an intrinsic accuracy of psi. Caller
// must establish forcing status and retain dataset identity. These autonomous
// defects on forced data are not controlled-predictor residuals.
struct EigenfunctionValidationResult {
    std::vector<EigenfunctionModeValidation> modes;
    std::size_t samples{};
    std::size_t stateDimension{}; // compatibility: dimension of supplied values
    std::size_t outsideTrainingSamples{};
    double outsideTrainingFraction{}, maximumRelativeExcursion{};
    bool trainingDomainAvailable{};
    EigenfunctionValidationInput inputKind{EigenfunctionValidationInput::physicalState};
    std::size_t inputDimension() const noexcept { return stateDimension; }
};
EigenfunctionValidationResult validateEigenfunctions(
    const KoopmanModel &model, const std::vector<double> &time, const std::vector<double> &values,
    std::size_t samples, std::size_t inputDimension, EigenfunctionValidationInput input,
    std::vector<double> horizons = {}, const std::vector<double> &trainingMinimum = {},
    const std::vector<double> &trainingMaximum = {});
EigenfunctionValidationResult
validateEigenfunctions(const PolynomialEigenfunctionSet &eigenfunctions,
                       const std::vector<double> &time, const std::vector<double> &coordinates,
                       std::vector<double> horizons = {},
                       const std::vector<double> &trainingMinimum = {},
                       const std::vector<double> &trainingMaximum = {});
} // namespace dynsyslab::core
