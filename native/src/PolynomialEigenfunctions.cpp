#include "PolynomialEvaluation.hpp"
#include <algorithm>
#include <cmath>
#include <dynsyslab/core/PolynomialEigenfunctions.hpp>
#include <limits>
#include <set>
#include <stdexcept>
namespace dynsyslab::core {
std::size_t PolynomialDictionary::featureCount() const noexcept {
    return coordinateDimension ? exponents.size() / coordinateDimension : 0;
}
void PolynomialDictionary::validate() const {
    if (!coordinateDimension || exponents.empty() || exponents.size() % coordinateDimension ||
        center.size() != coordinateDimension || scale.size() != coordinateDimension)
        throw std::invalid_argument("Polynomial dictionary dimensions are inconsistent");
    if (!std::all_of(center.begin(), center.end(), [](double v) { return std::isfinite(v); }) ||
        !std::all_of(scale.begin(), scale.end(),
                     [](double v) { return std::isfinite(v) && v > 0; }))
        throw std::invalid_argument(
            "Polynomial normalization needs finite centers and positive scales");
    std::set<std::vector<std::size_t>> unique;
    for (std::size_t f = 0; f < featureCount(); ++f) {
        const auto begin = exponents.begin() + f * coordinateDimension;
        if (!unique.emplace(begin, begin + coordinateDimension).second)
            throw std::invalid_argument("Duplicate polynomial monomial");
    }
}
struct PolynomialEigenfunctionSet::Data {
    PolynomialDictionary dictionary;
    std::vector<std::complex<double>> coefficients;
    EigenfunctionSpectrum spectrum;
    EigenfunctionCoordinates coordinates;
};
PolynomialEigenfunctionSet::PolynomialEigenfunctionSet(PolynomialDictionary d,
                                                       std::vector<std::complex<double>> c,
                                                       EigenfunctionSpectrum s,
                                                       EigenfunctionCoordinates q) {
    if (d.center.empty())
        d.center.assign(d.coordinateDimension, 0.);
    if (d.scale.empty())
        d.scale.assign(d.coordinateDimension, 1.);
    d.validate();
    q.validate();
    const auto features = d.featureCount();
    if (q.dimension != d.coordinateDimension ||
        s.size() > std::numeric_limits<std::size_t>::max() / features ||
        c.size() != features * s.size())
        throw std::invalid_argument("Polynomial coefficient/spectrum/coordinate dimensions differ");
    if (!std::all_of(c.begin(), c.end(),
                     [](auto z) { return std::isfinite(z.real()) && std::isfinite(z.imag()); }))
        throw std::invalid_argument("Polynomial coefficients must be finite");
    for (std::size_t j = 0; j < s.size(); ++j)
        if (std::all_of(c.begin() + j * features, c.begin() + (j + 1) * features,
                        [](auto z) { return z == std::complex<double>{}; }))
            throw std::invalid_argument(
                "Identically zero polynomial is not an eigenfunction candidate");
    data_ =
        std::make_shared<const Data>(Data{std::move(d), std::move(c), std::move(s), std::move(q)});
}
const PolynomialDictionary &PolynomialEigenfunctionSet::dictionary() const noexcept {
    return data_->dictionary;
}
const std::vector<std::complex<double>> &PolynomialEigenfunctionSet::coefficients() const noexcept {
    return data_->coefficients;
}
const EigenfunctionSpectrum &PolynomialEigenfunctionSet::spectrum() const noexcept {
    return data_->spectrum;
}
const EigenfunctionCoordinates &PolynomialEigenfunctionSet::coordinates() const noexcept {
    return data_->coordinates;
}
std::size_t PolynomialEigenfunctionSet::coordinateDimension() const noexcept {
    return dictionary().coordinateDimension;
}
std::size_t PolynomialEigenfunctionSet::eigenfunctionCount() const noexcept {
    return spectrum().size();
}
std::complex<double> PolynomialEigenfunctionSet::evaluate(std::size_t mode, const double *q) const {
    if (!q || mode >= eigenfunctionCount())
        throw std::invalid_argument("Invalid polynomial evaluation request");
    const auto &d = dictionary();
    return detail::polynomialValue(d.coordinateDimension, d.featureCount(), d.exponents, d.center,
                                   d.scale, coefficients().data() + mode * d.featureCount(), q);
}
std::vector<std::complex<double>> PolynomialEigenfunctionSet::evaluate(const double *q) const {
    std::vector<std::complex<double>> values;
    values.reserve(eigenfunctionCount());
    for (std::size_t j = 0; j < eigenfunctionCount(); ++j)
        values.push_back(evaluate(j, q));
    return values;
}
std::vector<std::complex<double>> PolynomialEigenfunctionSet::gradient(std::size_t mode,
                                                                       const double *q) const {
    if (!q || mode >= eigenfunctionCount())
        throw std::invalid_argument("Invalid polynomial gradient request");
    const auto &d = dictionary();
    return detail::polynomialGradient(d.coordinateDimension, d.featureCount(), d.exponents,
                                      d.center, d.scale,
                                      coefficients().data() + mode * d.featureCount(), q);
}
PolynomialEigenfunctionSet
PolynomialEigenfunctionSet::extract(const std::vector<std::size_t> &modes) const {
    auto selected = spectrum().extract(modes);
    std::vector<std::complex<double>> c;
    const auto n = dictionary().featureCount();
    for (auto j : modes)
        c.insert(c.end(), coefficients().begin() + j * n, coefficients().begin() + (j + 1) * n);
    return {dictionary(), std::move(c), std::move(selected), coordinates()};
}
PolynomialEigenfunctionSet
PolynomialEigenfunctionSet::withCoordinates(EigenfunctionCoordinates q) const {
    return {dictionary(), coefficients(), spectrum(), std::move(q)};
}
PolynomialEigenfunctionView PolynomialEigenfunctionSet::eigenfunction(std::size_t mode) const {
    if (mode >= eigenfunctionCount())
        throw std::out_of_range("Eigenfunction mode index");
    return {*this, mode};
}
} // namespace dynsyslab::core
