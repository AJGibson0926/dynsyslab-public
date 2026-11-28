#pragma once

#include <dynsyslab/core/Json.hpp>
#include <cstddef>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace dynsyslab::core {

// Immutable real expression DAG. Derivatives are expressions too, so repeated
// Lie derivatives remain exact formulas rather than nested finite differences. A node used several
// times (a named definition) is evaluated, differentiated and serialized once per call: JSON is a
// tree (analytic-state-expression/v1) when nothing is shared, else a node table (v2).
class AnalyticExpression {
public:
    struct Node;
    static AnalyticExpression constant(double value, std::size_t dimension);
    static AnalyticExpression state(std::size_t dimension, std::size_t index);
    static AnalyticExpression fromJson(const json::Value& value);
    // Infix text over named variables, e.g. "r^2*cos(theta)" (see AnalyticExpressionText.cpp).
    static AnalyticExpression parse(const std::string& text, const std::vector<std::string>& variables);
    // Named constants (e.g. a system's parameters) are substituted by value.
    static AnalyticExpression parse(const std::string& text, const std::vector<std::string>& variables,
                                    const std::map<std::string, double>& constants);
    [[nodiscard]] json::Value toJson() const;
    [[nodiscard]] std::size_t stateDimension() const noexcept;
    [[nodiscard]] double evaluate(const double* state) const;
    [[nodiscard]] AnalyticExpression derivative(std::size_t axis) const;
    [[nodiscard]] AnalyticExpression flowDerivative(
        const std::vector<AnalyticExpression>& vectorField, std::size_t order = 1) const;
    [[nodiscard]] AnalyticExpression power(double exponent) const;
    // exp, log, sin, cos, tanh, sech, ellipk (K(m)), ellipe (E(m)).
    [[nodiscard]] AnalyticExpression unary(const std::string& operation) const;
    // Two-argument functions: atan2(this, x), ellipf (F(this|m)), ellipeinc (E(this|m)); elliptic
    // integrals use the parameter m = k^2 (EllipticIntegrals.hpp).
    [[nodiscard]] AnalyticExpression function(const std::string& operation, const AnalyticExpression& second) const;
    // Functions of several arguments: ellippi(phi, n, m) = Pi(n; phi|m), the elliptic integral of the
    // third kind (DLMF 19.2.7; n < 1), and polyroot(k, c_0, ..., c_d), the k-th largest real root of
    // c_0 x^d + ... + c_d (PolynomialRoots.hpp; k a nonnegative integer constant). Both are exact
    // functions with exact derivatives (DLMF 19.4; implicit differentiation of the root).
    [[nodiscard]] static AnalyticExpression call(const std::string& operation, std::vector<AnalyticExpression> arguments);
    // A registered special function (ExpressionFunctions.cpp: Bessel, Airy, gamma family, error and
    // exponential integrals, Lambert W, hypergeometric, Jacobi amplitude, Weierstrass p, ...) with its
    // constant parameters (an order, hypergeometric parameters, invariants) and expression arguments.
    [[nodiscard]] static AnalyticExpression special(const std::string& name, std::vector<AnalyticExpression> parameters,
                                                    std::vector<AnalyticExpression> arguments);
    [[nodiscard]] static AnalyticExpression special(const std::string& name, std::vector<double> parameters,
                                                    std::vector<AnalyticExpression> arguments);
    // True when the expression depends on the state (not only on constants and summation variables).
    [[nodiscard]] bool dependsOnState() const noexcept;
    // Series. index(slot) is the summation variable of slot `slot` (0..63); sum(slot, from, to, term)
    // is sum_{k = from}^{to} term(k), `to` an integer or +infinity. An infinite sum adds terms until five
    // in a row are below 1e-17 of the partial sum (relative) and is NaN if that does not happen within
    // 1e6 terms: convergence is checked, never assumed. Derivatives are taken term by term, valid for
    // series whose term derivatives converge uniformly (as for the power and Fourier series here).
    [[nodiscard]] static AnalyticExpression index(std::size_t dimension, std::size_t slot);
    [[nodiscard]] static AnalyticExpression sum(std::size_t slot, double from, double to, const AnalyticExpression& term);
    AnalyticExpression operator+(const AnalyticExpression& other) const;
    AnalyticExpression operator-(const AnalyticExpression& other) const;
    AnalyticExpression operator*(const AnalyticExpression& other) const;
    AnalyticExpression operator/(const AnalyticExpression& other) const;
    AnalyticExpression operator-() const;
    AnalyticExpression operator+(double value) const;
    AnalyticExpression operator-(double value) const;
    AnalyticExpression operator*(double value) const;
    AnalyticExpression operator/(double value) const;
private:
    explicit AnalyticExpression(std::shared_ptr<const Node> node);
    double evaluate(const double* state, std::unordered_map<const Node*, double>& memo, std::vector<double>& indices) const;
    static std::shared_ptr<const Node> make(Node node);  // with the summation variables it depends on
    static AnalyticExpression rebuild(const json::Value& node, std::size_t dimension, std::vector<AnalyticExpression> arguments);
    AnalyticExpression derivative(std::size_t axis, std::unordered_map<const Node*, AnalyticExpression>& memo) const;
    AnalyticExpression binary(const std::string&, const AnalyticExpression&) const;
    std::shared_ptr<const Node> node_;
};
AnalyticExpression operator+(double value, const AnalyticExpression& expression);
AnalyticExpression operator-(double value, const AnalyticExpression& expression);
AnalyticExpression operator*(double value, const AnalyticExpression& expression);
AnalyticExpression operator/(double value, const AnalyticExpression& expression);

}  // namespace dynsyslab::core
