#pragma once
// Truncated multivariate power series with complex coefficients, and the evaluation of an analytic
// state expression on them: the Taylor coefficients of f(x0 + W(p)) for the parameterization method
// (SpectralSubmanifold.hpp). Elementary functions compose through their Taylor series about the
// constant term (automatic differentiation by power-series arithmetic; Griewank and Walther,
// Evaluating Derivatives, 2nd ed. 2008, chapter 13; Haro et al., The Parameterization Method for
// Invariant Manifolds, Springer 2016, section 2.2).
#include <dynsyslab/core/Json.hpp>

#include <complex>
#include <cstddef>
#include <memory>
#include <vector>

namespace dynsyslab::core::detail {

using Complex = std::complex<double>;

// Monomials p^m in `variables` variables of total degree <= order, graded (degree 0 first, then
// lexicographically descending exponents within a degree), with their product table.
class MonomialTable {
public:
    MonomialTable(std::size_t variables, std::size_t order);
    [[nodiscard]] std::size_t size() const noexcept { return degrees_.size(); }
    [[nodiscard]] std::size_t variables() const noexcept { return variables_; }
    [[nodiscard]] std::size_t order() const noexcept { return order_; }
    [[nodiscard]] std::size_t degree(std::size_t i) const { return degrees_[i]; }
    [[nodiscard]] const std::size_t* exponents(std::size_t i) const { return exponents_.data() + i * variables_; }
    // Index of the monomial with these exponents; size() when its degree exceeds the order.
    [[nodiscard]] std::size_t index(const std::size_t* exponents) const;
    // Index of p^a p^b, or size() when it exceeds the order.
    [[nodiscard]] std::size_t product(std::size_t a, std::size_t b) const { return products_[a * size() + b]; }
    // First and one-past-last index of the monomials of degree k.
    [[nodiscard]] std::size_t begin(std::size_t k) const { return starts_[k]; }
    [[nodiscard]] std::size_t end(std::size_t k) const { return starts_[k + 1]; }

private:
    std::size_t variables_, order_;
    std::vector<std::size_t> exponents_, degrees_, starts_, products_;
};

class Series {
public:
    explicit Series(std::shared_ptr<const MonomialTable> table, Complex constant = 0);
    [[nodiscard]] static Series variable(std::shared_ptr<const MonomialTable> table, std::size_t index, Complex at = 0);
    [[nodiscard]] const MonomialTable& table() const { return *table_; }
    [[nodiscard]] const std::shared_ptr<const MonomialTable>& shared() const { return table_; }
    std::vector<Complex> c;  // one coefficient per monomial of the table

    Series operator+(const Series& b) const;
    Series operator-(const Series& b) const;
    Series operator*(const Series& b) const;
    Series operator*(Complex s) const;
    Series operator-() const;
    // The part of degree < k (k = 0 gives zero).
    [[nodiscard]] Series truncated(std::size_t k) const;

private:
    std::shared_ptr<const MonomialTable> table_;
};

[[nodiscard]] Series exp(const Series& u);
[[nodiscard]] Series log(const Series& u);
[[nodiscard]] Series sin(const Series& u);
[[nodiscard]] Series cos(const Series& u);
[[nodiscard]] Series power(const Series& u, double exponent);
[[nodiscard]] Series reciprocal(const Series& u);

// The series of an analytic-state-expression/v1 tree (AnalyticExpression::toJson) with state
// coordinate i replaced by `state[i]`. Throws for functions without a power series here (atan2,
// elliptic integrals).
[[nodiscard]] Series evaluate(const json::Value& expression, const std::vector<Series>& state);

}  // namespace dynsyslab::core::detail
