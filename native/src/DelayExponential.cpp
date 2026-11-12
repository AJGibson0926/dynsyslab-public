#include <boost/math/special_functions/gamma.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <cmath>
#include <dynsyslab/core/Delay.hpp>
#include <stdexcept>
namespace dynsyslab::core {
namespace {
using Real = boost::multiprecision::cpp_dec_float_50;
// Integral_0^t s^k exp(a*s) ds / k!, for a>0. All series terms are nonnegative.
Real growingIntegral(unsigned k, const Real &a, const Real &t, const Real &factorial) {
    if (t == 0)
        return 0;
    Real term = pow(t, k + 1) / (factorial * (k + 1)), sum = term;
    for (unsigned j = 0; j < 100000; ++j) {
        term *= a * t * Real(k + j + 1) / (Real(j + 1) * Real(k + j + 2));
        sum += term;
        if (abs(term) <= abs(sum) * Real("1e-45"))
            return sum;
    }
    throw std::domain_error("Delayed exponential series did not converge at this scale");
}
} // namespace
double linearDelayConstantHistory(double beta, double gamma, double delay, double history,
                                  double elapsed) {
    if (!std::isfinite(beta) || !std::isfinite(gamma) || !std::isfinite(delay) || delay <= 0 ||
        !std::isfinite(history) || !std::isfinite(elapsed))
        throw std::invalid_argument(
            "Delayed exponential requires finite arguments and positive delay");
    if (elapsed <= 0 || history == 0 || beta == gamma)
        return history;
    if (beta == 0) {
        const double value = history * std::exp(-gamma * elapsed);
        if (!std::isfinite(value))
            throw std::domain_error("Delayed exponential left the finite domain");
        return value;
    }
    if (elapsed / delay > 100000)
        throw std::domain_error("Delayed exponential exceeds 100000 exact method-of-steps terms");
    const unsigned steps = static_cast<unsigned>(std::floor(elapsed / delay));
    const Real b(beta), g(gamma), tau(delay), t(elapsed);
    Real power = 1, factorial = 1, gammaPower = g, result = 0, absoluteSum = 0;
    for (unsigned k = 0; k <= steps; ++k) {
        Real upper = t - k * tau;
        if (upper < 0)
            upper = 0; // only a rounded interval endpoint can cross zero
        Real lower = upper - tau;
        if (lower < 0)
            lower = 0;
        Real fundamental = power * pow(upper, k) * exp(-g * upper) / factorial;
        Real integral;
        if (gamma > 0) {
            Real a = g * lower, z = g * upper;
            // Complementary tails avoid subtraction of two rounded values near one.
            Real difference =
                a > k + 1
                    ? boost::math::gamma_q(Real(k + 1), a) - boost::math::gamma_q(Real(k + 1), z)
                    : boost::math::gamma_p(Real(k + 1), z) - boost::math::gamma_p(Real(k + 1), a);
            integral = difference / gammaPower;
        } else if (gamma == 0) {
            integral = (pow(upper, k + 1) - pow(lower, k + 1)) / (factorial * (k + 1));
        } else {
            integral =
                growingIntegral(k, -g, upper, factorial) - growingIntegral(k, -g, lower, factorial);
        }
        const Real contribution = fundamental + b * power * integral;
        result += contribution;
        absoluteSum += abs(fundamental) + abs(b * power * integral);
        power *= b;
        factorial *= k + 1;
        gammaPower *= g;
    }
    if (absoluteSum > (abs(result) > 1 ? abs(result) : Real(1)) * Real("1e30"))
        throw std::domain_error("Delayed exponential cancellation exceeds its precision budget");
    const double value = (Real(history) * result).convert_to<double>();
    if (!std::isfinite(value))
        throw std::domain_error("Delayed exponential left the finite domain");
    return value;
}
} // namespace dynsyslab::core
