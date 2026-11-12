#pragma once
#include <complex>
#include <vector>
namespace dynsyslab::core::detail {
// Shared scalar evaluation kernel; same operation ordering as the historical
// SINDy evaluator. Fitting and runtime representations do not duplicate it.
inline double integerPower(double x, std::size_t exponent) {
    double result = 1.;
    for (std::size_t i = 0; i < exponent; ++i)
        result *= x;
    return result;
}
inline std::complex<double> polynomialValue(std::size_t dimension, std::size_t features,
                                            const std::vector<std::size_t> &exponents,
                                            const std::vector<double> &center,
                                            const std::vector<double> &scale,
                                            const std::complex<double> *coefficients,
                                            const double *q) {
    std::complex<double> result{};
    for (std::size_t f = 0; f < features; ++f) {
        double monomial = 1.;
        for (std::size_t a = 0; a < dimension; ++a)
            monomial *= integerPower((q[a] - center[a]) / scale[a], exponents[f * dimension + a]);
        result += coefficients[f] * monomial;
    }
    return result;
}
inline std::vector<std::complex<double>>
polynomialGradient(std::size_t dimension, std::size_t features,
                   const std::vector<std::size_t> &exponents, const std::vector<double> &center,
                   const std::vector<double> &scale, const std::complex<double> *coefficients,
                   const double *q) {
    std::vector<std::complex<double>> result(dimension);
    for (std::size_t f = 0; f < features; ++f)
        for (std::size_t d = 0; d < dimension; ++d) {
            const auto count = exponents[f * dimension + d];
            if (!count)
                continue;
            double derivative = double(count) / scale[d];
            for (std::size_t a = 0; a < dimension; ++a) {
                auto e = exponents[f * dimension + a];
                if (a == d)
                    --e;
                derivative *= integerPower((q[a] - center[a]) / scale[a], e);
            }
            result[d] += coefficients[f] * derivative;
        }
    return result;
}
} // namespace dynsyslab::core::detail
