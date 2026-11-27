#include <dynsyslab/core/EigenfunctionRefit.hpp>

#include <dynsyslab/numerics/MatrixFunctions.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::core {

std::vector<std::complex<double>> refitAtMultiplier(const std::vector<double>& now, const std::vector<double>& next,
                                                    std::size_t pairs, std::size_t terms, std::complex<double> mu,
                                                    double* residual) {
    if (now.size() != pairs * terms || next.size() != pairs * terms || pairs < terms)
        throw std::invalid_argument("refitAtMultiplier needs at least as many sample pairs as dictionary terms");
    // D^H D for D = next - mu now (complex, terms x terms); its smallest eigenvector is the fit.
    std::vector<std::complex<double>> gram(terms * terms, 0.0);
    std::vector<std::complex<double>> row(terms);
    for (std::size_t p = 0; p < pairs; ++p) {
        for (std::size_t k = 0; k < terms; ++k) row[k] = next[p * terms + k] - mu * now[p * terms + k];
        for (std::size_t i = 0; i < terms; ++i)
            for (std::size_t j = 0; j < terms; ++j) gram[i * terms + j] += std::conj(row[i]) * row[j];
    }
    const auto eigen = numerics::hermitianEigensystem(gram, terms);
    std::vector<std::complex<double>> coefficients(terms);
    for (std::size_t k = 0; k < terms; ++k) coefficients[k] = eigen.eigenvectors[k * terms];  // column 0
    if (residual) *residual = std::sqrt(std::max(0.0, eigen.eigenvalues.front()) / double(pairs));
    return coefficients;
}

std::complex<double> gaugeFactor(const std::vector<std::complex<double>>& values,
                                 const std::vector<std::complex<double>>& reference) {
    if (values.size() != reference.size() || values.empty()) throw std::invalid_argument("gauge needs aligned samples");
    std::complex<double> numerator = 0.0;
    double denominator = 0.0;
    for (std::size_t k = 0; k < values.size(); ++k) {
        numerator += std::conj(values[k]) * reference[k];
        denominator += std::norm(values[k]);
    }
    if (!(denominator > 0.0)) throw std::invalid_argument("gauge: the eigenfunction vanishes on the samples");
    return numerator / denominator;
}

std::size_t pruneCoefficients(std::vector<std::complex<double>>& coefficients, double threshold) {
    double largest = 0.0;
    for (const auto& c : coefficients) largest = std::max(largest, std::abs(c));
    std::size_t removed = 0;
    for (auto& c : coefficients)
        if (std::abs(c) < threshold * largest && c != 0.0) { c = 0.0; ++removed; }
    return removed;
}

}  // namespace dynsyslab::core
