#include <dynsyslab/core/EigenfunctionFit.hpp>
namespace dynsyslab::core {
EigenfunctionFitResult fitSindyEigenfunctions(const std::vector<double> &values,
                                              std::size_t samples, const std::vector<double> &time,
                                              const std::vector<std::complex<double>> &rates,
                                              FunctionDomain coordinates,
                                              const std::vector<double> &derivatives,
                                              const KoopmanSINDyOptions &options) {
    coordinates.validate();
    auto legacy = fitKoopmanEigenfunctions(values, samples, coordinates.dimension, time, rates,
                                           derivatives, options);
    auto diagnostics = fitDiagnosticsFromLegacy(legacy);
    diagnostics.algorithm = "sindy-generator";
    diagnostics.residualDefinition = "Dpsi(q) qdot - lambda psi(q)";
    auto fit = polynomialEigenfunctionsFromLegacy(legacy);
    return {fit.functions.withCoordinates(std::move(coordinates)), std::move(fit.spectrum), std::move(diagnostics)};
}
} // namespace dynsyslab::core
