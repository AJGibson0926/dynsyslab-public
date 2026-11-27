#pragma once
// Held-out quality of a Koopman readout y ~= C psi: fit on a training split, score on the held-out
// split (or on an independent validation set) by R^2 and RMSE normalised by the output's
// standard deviation per output, with the condition number of the training design (Psi^H Psi)
// (cross-validation: Stone 1974, "Cross-validatory choice and assessment of statistical
// predictions", J. R. Stat. Soc. B 36:111-147). A gate compares the scores with thresholds.
#include <dynsyslab/core/Reconstruction.hpp>

#include <complex>
#include <cstddef>
#include <string>
#include <vector>

namespace dynsyslab::core {

struct ReconstructionQualityPolicy {
    double validationFraction{0.25};
    std::string split{"interleaved"};  // interleaved | tail
    double minimumR2{0.99};
    double maximumNormalizedRmse{0.10};
    double maximumConditionNumber{1e12};
};

struct ReconstructionQuality {
    ReconstructionResult fit;           // fitted on the training rows only
    std::vector<double> r2;             // per output, on the held-out rows
    std::vector<double> normalizedRmse; // per output
    double conditionNumber{};
    std::size_t trainingSamples{}, validationSamples{};
    bool passed{};
    std::vector<std::string> failures;
};

// psi: samples x modes, y: samples x outputs (row-major). With validation arrays given, all rows
// train and the validation rows score.
[[nodiscard]] ReconstructionQuality assessReconstruction(
    const std::vector<std::complex<double>>& psi, const std::vector<std::complex<double>>& y,
    std::size_t samples, std::size_t modes, std::size_t outputs, double ridge,
    const ReconstructionQualityPolicy& policy,
    const std::vector<std::complex<double>>& validationPsi = {},
    const std::vector<std::complex<double>>& validationY = {}, std::size_t validationSamples = 0);

}  // namespace dynsyslab::core
