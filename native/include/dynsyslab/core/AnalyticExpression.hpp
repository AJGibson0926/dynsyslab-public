#pragma once

#include <dynsyslab/core/Json.hpp>
#include <cstddef>
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
    [[nodiscard]] json::Value toJson() const;
    [[nodiscard]] std::size_t stateDimension() const noexcept;
    [[nodiscard]] double evaluate(const double* state) const;
    [[nodiscard]] AnalyticExpression derivative(std::size_t axis) const;
    [[nodiscard]] AnalyticExpression flowDerivative(
        const std::vector<AnalyticExpression>& vectorField, std::size_t order = 1) const;
    [[nodiscard]] AnalyticExpression power(double exponent) const;
    [[nodiscard]] AnalyticExpression unary(const std::string& operation) const;
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
