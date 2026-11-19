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
// Classical observability matrix of a K-ROM in discrete modal coordinates:
//   O = [C; C M; C M^2; ...; C M^(blocks-1)],   M = diag(mu_j), mu_j = exp(lambda_j dt),
// with C the complex readout (outputs x modes) of `measurement`, optionally restricted to
// `outputRows`. One coordinate per eigenfunction, so a conjugate pair counts as two complex
// modes. blocks = 0 uses the mode count; dt = 0 uses the K-ROM's sample interval.
// Singular values (descending) and right singular vectors (observability directions,
// modes x k, column k) come from a direct complex SVD of O.
struct ModalObservabilityMatrix {
    std::vector<std::complex<double>> matrix;      // (blocks * outputs) x modes, row-major
    std::vector<double> singularValues;
    std::vector<std::complex<double>> directions;  // modes x rank, column k = direction k
    std::vector<std::complex<double>> multipliers; // mu_j
    std::size_t modes{}, outputs{}, blocks{};
    double sampleTime{};
};
[[nodiscard]] ModalObservabilityMatrix modalObservabilityMatrix(
    const KoopmanModel& model, const std::string& measurement, double sampleTime = 0.0,
    std::size_t blocks = 0, const std::vector<std::size_t>& outputRows = {});
} // namespace dynsyslab::core
