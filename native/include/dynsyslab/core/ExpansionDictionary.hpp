#pragma once
// The dictionary of an eigenfunction expansion psi_j(q) = sum_k c_jk phi_k(q): an ordered collection
// of basis functions phi_k of any kind. Methods (SINDy, EDMD, EDMDc, DMD, spectral fits, or a person
// writing coefficients by hand) differ only in which terms they choose and how they supply c.
//
// Term kinds, all evaluated on q' = (q - center) / scale (identity when center/scale are empty):
//   coordinate                q'_index, one coordinate (a linear term stored by its index alone,
//                             so long linear maps, such as a window map, stay small)
//   monomial                  prod_i q'_i^n_i
//   chebyshev, laguerre,      prod_i P_{n_i}(q'_i) (NIST DLMF ch. 18); fourier_hermite: Fourier in
//   hermite, fourier_hermite  q'_0 times Hermite functions in q'_1
//   gaussian, laplacian,      radial kernels of |q' - node'| with width, smoothness, amplitude
//   lorentzian, matern
//   sine, cosine              sin / cos(w . q' + phase)
//   exponential               exp(w . q')
//   partition                 1 if `node` is the nearest of `nodes`, else 0 (Voronoi indicator)
//   expression                any analytic expression of the raw q
#include <dynsyslab/core/AnalyticExpression.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/SpectralKROM.hpp>

#include <complex>
#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace dynsyslab::core {

class Observable;

// A table of monomials sharing one normalization (the SINDy / EDMD polynomial library).
struct PolynomialDictionary {
    std::size_t coordinateDimension{};
    std::vector<std::size_t> exponents; // feature x coordinate
    std::vector<double> center, scale;
    std::size_t featureCount() const noexcept;
    void validate() const;
    // Reuse this same basis as a real candidate library for dynamics or other regressions.
    Observable asObservable() const;
};

struct DictionaryTerm {
    std::string kind;
    std::vector<double> center, scale;
    std::vector<std::size_t> orders;   // product kinds: one order per coordinate
    std::vector<double> vector;        // radial node, trig frequency or exponential rate (raw units for nodes)
    double phase{0}, width{1}, smoothness{2}, amplitude{1};
    std::vector<double> nodes;         // partition: all nodes (raw), this term is node `index`
    std::size_t index{};
    std::shared_ptr<const AnalyticExpression> expression;
};

class ExpansionDictionary {
public:
    ExpansionDictionary(std::size_t dimension, std::vector<DictionaryTerm> terms);
    [[nodiscard]] static ExpansionDictionary monomials(PolynomialDictionary table);
    [[nodiscard]] static ExpansionDictionary fromSpectral(const SpectralDictionary& dictionary);

    [[nodiscard]] std::size_t dimension() const noexcept { return dimension_; }
    [[nodiscard]] std::size_t size() const noexcept { return terms_.size(); }
    [[nodiscard]] const std::vector<DictionaryTerm>& terms() const noexcept { return terms_; }
    // Distinct term kinds, in first-use order.
    [[nodiscard]] std::vector<std::string> kinds() const;
    // Non-null when every term is a monomial sharing one normalization.
    [[nodiscard]] const PolynomialDictionary* monomialTable() const noexcept { return table_ ? &*table_ : nullptr; }

    [[nodiscard]] std::vector<double> evaluate(const double* q) const;   // size()
    [[nodiscard]] std::vector<double> gradient(const double* q) const;   // size() x dimension()
    [[nodiscard]] bool exactGradient() const noexcept;
    [[nodiscard]] std::complex<double> combine(const std::complex<double>* c, const double* q) const;
    [[nodiscard]] std::vector<std::complex<double>> combineGradient(const std::complex<double>* c, const double* q) const;

    [[nodiscard]] json::Value describe() const;
    [[nodiscard]] static ExpansionDictionary fromJson(const json::Value& value);

private:
    std::size_t dimension_{};
    std::vector<DictionaryTerm> terms_;
    std::optional<PolynomialDictionary> table_;
};

}  // namespace dynsyslab::core
