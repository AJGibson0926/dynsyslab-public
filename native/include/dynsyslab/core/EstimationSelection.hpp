#pragma once
#include <dynsyslab/core/Estimation.hpp>
namespace dynsyslab::core {
struct EstimationSelectionOptions {
    std::string measurement{"echo"}, output{"state"};
    std::vector<double> outputScales; // physical accuracy units, never silently mixed
    double measurementStandardDeviation{0.01}; // independent sample noise
    double horizon{1.0}, sampleTime{0.01};
    double dynamicsWeight{1.0}, measurementMismatchWeight{1.0};
};
struct EstimationSubsetAssessment {
    double stateReconstructionMse{}, measurementReconstructionMse{};
    double posteriorStateNoiseMse{}, autonomousPropagationMse{}, score{};
    std::vector<double> modeScales;
    ObservabilityResult observability;
};
// Assess a proposed subset using its retained C columns (no reconstruction
// refitting). Psi and reference outputs are aligned sample rows. Noise score
// assumes a unit isotropic prior AFTER modal RMS normalization; it is a
// finite-window linear Gaussian benchmark, not a guarantee of KLQE accuracy.
[[nodiscard]] EstimationSubsetAssessment assessEstimationSubset(
    const KoopmanModel& subset,const std::vector<double>& time,
    const std::vector<std::complex<double>>& referencePsi,
    const std::vector<double>& actualOutput,const std::vector<double>& actualMeasurement,
    const EstimationSelectionOptions& options);
} // namespace dynsyslab::core
