#pragma once

#include <dynsyslab/core/Json.hpp>
#include <cstddef>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace dynsyslab::core {

// Immutable real expression DAG. Derivatives are expressions too, so repeated
// Lie derivatives remain exact formulas rather than nested finite differences.
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
    AnalyticExpression binary(const std::string&, const AnalyticExpression&) const;
    std::shared_ptr<const Node> node_;
};
AnalyticExpression operator+(double value, const AnalyticExpression& expression);
AnalyticExpression operator-(double value, const AnalyticExpression& expression);
AnalyticExpression operator*(double value, const AnalyticExpression& expression);
AnalyticExpression operator/(double value, const AnalyticExpression& expression);

}  // namespace dynsyslab::core
