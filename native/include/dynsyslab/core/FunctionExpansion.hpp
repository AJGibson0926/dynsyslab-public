#pragma once
// The expansion representation: g_r(q) = sum_k A_rk phi_k(q), r = 0..outputs-1, in the one
// universal dictionary (ExpansionDictionary) on a domain q. It holds dictionary terms and
// coefficients only. Observables and K-ROMs point at it (Observable::representation()); an
// eigenvalue belongs to the Eigenfunction, never to the expansion. Methods (SINDy, EDMD, DMD,
// spectral fits, library entries, jobs) differ only in which terms and coefficients they supply.
#include <dynsyslab/core/ExpansionDictionary.hpp>
#include <dynsyslab/core/FunctionRepresentation.hpp>

#include <complex>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace dynsyslab::core {

class FunctionExpansion final : public FunctionRepresentation {
public:
    FunctionExpansion(std::shared_ptr<const ExpansionDictionary> dictionary, std::vector<std::complex<double>> coefficients,
                        FunctionDomain domain, std::vector<std::string> outputNames = {});
    // A dictionary of monomials sharing one normalization.
    FunctionExpansion(PolynomialDictionary monomials, std::vector<std::complex<double>> coefficients, FunctionDomain domain);

    // Real analytic expressions over the domain coordinates, decomposed into dictionary terms
    // (polynomial parts become monomials; other parts become expression terms).
    [[nodiscard]] static FunctionExpansion fromExpressions(const std::vector<AnalyticExpression>& expressions,
                                                             FunctionDomain domain, std::vector<std::string> outputNames = {});

    const ExpansionDictionary& dictionary() const noexcept;
    std::shared_ptr<const ExpansionDictionary> dictionaryHandle() const noexcept;
    // The monomial table when every term is a monomial with one normalization; throws otherwise.
    const PolynomialDictionary& monomialTable() const;
    bool monomial() const noexcept;
    const std::vector<std::complex<double>>& coefficients() const noexcept;  // outputs x terms
    std::string kind() const override { return "expansion"; }
    const FunctionDomain& coordinates() const noexcept override;
    const std::vector<std::string>& outputNames() const noexcept;
    std::size_t outputCount() const noexcept override;
    bool exactJacobian() const noexcept override { return dictionary().exactGradient(); }

    std::complex<double> evaluate(std::size_t output, const double* coordinates) const;
    std::vector<std::complex<double>> evaluate(const double* coordinates) const override;
    std::vector<std::complex<double>> gradient(std::size_t output, const double* coordinates) const;
    std::vector<std::complex<double>> jacobian(const double* coordinates) const override;  // outputs x dimension

    // The same functions written in a larger dictionary that contains every term of this one.
    FunctionExpansion inDictionary(std::shared_ptr<const ExpansionDictionary> larger) const;
    FunctionExpansion select(const std::vector<std::size_t>& outputs) const;
    // The same functions on a relabelled domain (same dimension), e.g. with its state map resolved.
    FunctionExpansion withCoordinates(FunctionDomain coordinates) const;
    std::shared_ptr<const FunctionRepresentation> selectOutputs(const std::vector<std::size_t>& outputs) const override;
    std::shared_ptr<const FunctionRepresentation> onDomain(FunctionDomain coordinates) const override {
        return std::make_shared<const FunctionExpansion>(withCoordinates(std::move(coordinates)));
    }
    // Normalized monomials rewritten as raw monomials (binomial expansion), so that equal functions
    // share terms; other terms unchanged. Projections work on canonical forms.
    FunctionExpansion canonical() const;
    // Closed forms of the outputs as real expressions (real part, imaginary part per output), when
    // every term with a nonzero coefficient has one (all kinds except Matern and partitions).
    std::optional<std::vector<AnalyticExpression>> realExpressions() const;

private:
    struct Data;
    std::shared_ptr<const Data> data_;
};

// A dictionary term in closed form over the domain coordinates, if it has one.
[[nodiscard]] std::optional<AnalyticExpression> dictionaryTermExpression(const DictionaryTerm& term, std::size_t dimension);

// The union of two dictionaries on one domain: `first`'s terms in order, then `second`'s new terms.
// Terms are identified by value (kind, normalization and parameters).
[[nodiscard]] std::shared_ptr<const ExpansionDictionary> mergeDictionaries(const ExpansionDictionary& first,
                                                                           const ExpansionDictionary& second);

}  // namespace dynsyslab::core
