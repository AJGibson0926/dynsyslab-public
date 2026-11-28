#pragma once
// The expression language's special functions (ExpressionFunctions.cpp): each a value and the exact
// partial derivative with respect to each argument, itself an expression. Leading "parameters" (a
// Bessel order, hypergeometric parameters, Weierstrass invariants) are expressions that do not depend on
// the state, though they may on a summation variable (a series over the order); the remaining arguments
// are any expressions. A partial that is not available throws when that argument varies.
#include <dynsyslab/core/AnalyticExpression.hpp>

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

namespace dynsyslab::core::detail {

struct SpecialFunction {
    std::size_t parameters{}, arguments{};
    std::function<double(const std::vector<double>& parameters, const std::vector<double>& arguments)> value;
    // d f / d argument i, given the parameters, the arguments and f itself.
    std::function<AnalyticExpression(const std::vector<AnalyticExpression>& parameters, const std::vector<AnalyticExpression>& arguments,
                                     const AnalyticExpression& self, std::size_t i)> partial;
};

// The registered function `name`, or null.
[[nodiscard]] const SpecialFunction* specialFunction(const std::string& name);
// Names of every registered function, for documentation and error messages.
[[nodiscard]] std::vector<std::string> specialFunctionNames();
// Expressions the parser builds from simpler ones (sn, cn, dn, the orthogonal polynomials, the
// theta functions as series, using index slot `slot`): false when `name` is not one.
[[nodiscard]] bool expandSugar(const std::string& name, const std::vector<AnalyticExpression>& arguments, std::size_t slot,
                               AnalyticExpression& result);

}  // namespace dynsyslab::core::detail
