#include <dynsyslab/core/SINDy.hpp>
namespace dynsyslab::core {
EigenfunctionFitResult fitSindyEigenfunctions(const std::vector<double> &values,
                                              std::size_t samples, const std::vector<double> &time,
                                              const std::vector<std::complex<double>> &rates,
                                              FunctionDomain coordinates,
                                              const std::vector<double> &derivatives,
                                              const KoopmanSINDyOptions &options) {
    coordinates.validate();
    auto fit = fitKoopmanEigenfunctions(values, samples, coordinates.dimension, time, rates, derivatives, options);
    fit.functions = fit.functions.withCoordinates(std::move(coordinates));
    return fit;
}
} // namespace dynsyslab::core
