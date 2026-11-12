#pragma once
#include <dynsyslab/core/Observability.hpp>
#include <complex>
namespace dynsyslab::core {
struct ModalObservabilityResult {
    std::vector<std::complex<double>> gramian; // mode x mode, Hermitian
    std::vector<std::complex<double>> directions; // columns, strongest first
    std::vector<double> gramianEigenvalues, singularValues, coordinateScales;
    std::size_t modes{}, rank{};
    double rankTolerance{};
};
// Closed-form diagonal-dynamics Gramian: O(modes^2) assembly. Scales specify
// psi_j = scale_j * normalized_psi_j; output units/noise whitening are the
// caller's responsibility. This complex-coordinate analysis does not silently
// identify conjugate pairs as independent physical degrees of freedom.
[[nodiscard]] ModalObservabilityResult analyzeModalObservability(
    const KoopmanModel& model,const std::string& measurement,
    const ObservabilityOptions& options,
    const std::vector<double>& coordinateScales = {});
} // namespace dynsyslab::core
