#include <dynsyslab/core/SINDy.hpp>
#include "PolynomialEvaluation.hpp"
#include <stdexcept>

namespace dynsyslab::core {
namespace {
void monomials(std::size_t dimension, std::size_t remaining, std::size_t start,
               std::vector<std::size_t>& exponent, std::vector<std::size_t>& output) {
    if (!remaining) { output.insert(output.end(), exponent.begin(), exponent.end()); return; }
    for (std::size_t index = start; index < dimension; ++index) {
        ++exponent[index]; monomials(dimension, remaining - 1, index, exponent, output);
        --exponent[index];
    }
}

}
std::vector<std::size_t> polynomialExponentTable(
    std::size_t dimension, std::size_t degree, bool includeConstant) {
    if (!dimension || !degree)
        throw std::invalid_argument(
            "polynomial exponent table requires positive dimension and degree");
    std::vector<std::size_t> result, exponent(dimension, 0);
    if (includeConstant) result.insert(result.end(), dimension, 0);
    for (std::size_t order = 1; order <= degree; ++order)
        monomials(dimension, order, 0, exponent, result);
    return result;
}


std::complex<double> evaluatePolynomialEigenfunction(const KoopmanSINDyResult& model,
                                                     std::size_t eigenfunction,
                                                     const double* state) {
    if (eigenfunction >= model.continuousEigenvalues.size()) throw std::out_of_range("eigenfunction index");
    return detail::polynomialValue(model.stateDimension,model.featureCount,model.exponents,
        model.center,model.scale,model.coefficients.data()+eigenfunction*model.featureCount,state);
}

std::vector<std::complex<double>> gradientPolynomialEigenfunction(
    const KoopmanSINDyResult& model, std::size_t eigenfunction, const double* state) {
    if (eigenfunction >= model.continuousEigenvalues.size()) throw std::out_of_range("eigenfunction index");
    return detail::polynomialGradient(model.stateDimension,model.featureCount,model.exponents,
        model.center,model.scale,model.coefficients.data()+eigenfunction*model.featureCount,state);
}

} // namespace dynsyslab::core
