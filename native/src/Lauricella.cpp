#include <boost/math/quadrature/tanh_sinh.hpp>
#include <cmath>
#include <dynsyslab/core/SpecialFunctions.hpp>
#include <stdexcept>
namespace dynsyslab::core::special {
double lauricellaFD(double a, const std::vector<double> &b, double c,
                    const std::vector<double> &x) {
    if (!std::isfinite(a) || !std::isfinite(c) || !(c > a && a > 0) || b.size() != x.size())
        throw std::invalid_argument(
            "Real Lauricella integral requires c>a>0 and matching b/x vectors");
    double endpoint = c - a;
    for (std::size_t i = 0; i < x.size(); ++i) {
        if (!std::isfinite(b[i]) || !std::isfinite(x[i]) || x[i] > 1)
            throw std::invalid_argument("Real Lauricella arguments must be finite and x<=1");
        if (x[i] == 1)
            endpoint -= b[i];
    }
    if (!(endpoint > 0))
        throw std::domain_error("Divergent Lauricella endpoint");
    if (x.empty())
        return 1;
    const double normalization = std::lgamma(c) - std::lgamma(a) - std::lgamma(c - a);
    if (!std::isfinite(normalization))
        throw std::domain_error("Lauricella normalization overflow");
    boost::math::quadrature::tanh_sinh<double> quadrature(12);
    double error = 0, l1 = 0;
    const double value = quadrature.integrate(
        [&](double u, double complement) {
            // Boost supplies the signed distance to the nearest endpoint, avoiding
            // loss of 1-u at the singular right endpoint.
            const double v = complement > 0 ? complement : 1 + complement;
            if (complement < 0)
                u = -complement;
            double logarithm = normalization + (a - 1) * std::log(u) + (c - a - 1) * std::log(v);
            for (std::size_t i = 0; i < x.size(); ++i) {
                const double factor = x[i] == 1 ? std::log(v)
                                      : u > .5  ? std::log((1 - x[i]) + x[i] * v)
                                                : std::log1p(-x[i] * u);
                logarithm -= b[i] * factor;
            }
            return std::exp(logarithm);
        },
        0.0, 1.0, 2e-12, &error, &l1);
    if (!std::isfinite(value) || !std::isfinite(error) || error > 2e-10 * std::max(1.0, l1))
        throw std::runtime_error("Lauricella integral did not converge in the finite real domain");
    return value;
}
} // namespace dynsyslab::core::special
