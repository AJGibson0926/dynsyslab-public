#pragma once
// Real roots of a real polynomial, as an algebraic function of its coefficients: the expression
// language's polyroot (AnalyticExpression.hpp), whose gradient follows by implicit differentiation,
// dx/dc_i = -x^(d-i) / P'(x). All roots are found together by the Aberth-Ehrlich iteration (Aberth,
// Math. Comp. 27, 1973) in long double (extended where the platform has it) and the real ones polished by Newton's method on the
// real polynomial, so a simple real root is accurate to rounding. A root is real when its imaginary
// part is within a relative 1e-8 of its magnitude; near-double roots are therefore ill-conditioned,
// as they are for any method.
#include <cstddef>
#include <vector>

namespace dynsyslab::core {

// The real roots of c[0] x^d + c[1] x^(d-1) + ... + c[d], largest first (c[0] != 0, d >= 1).
[[nodiscard]] std::vector<double> realPolynomialRoots(const std::vector<double>& coefficients);
// The k-th largest real root (k = 0 the largest); NaN when there are fewer than k + 1.
[[nodiscard]] double realPolynomialRoot(std::size_t k, const std::vector<double>& coefficients);

}  // namespace dynsyslab::core
