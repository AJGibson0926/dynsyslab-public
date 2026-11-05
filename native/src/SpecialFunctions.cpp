#include <algorithm>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/jacobi_elliptic.hpp>
#include <boost/math/special_functions/jacobi_theta.hpp>
#include <cmath>
#include <dynsyslab/core/SpecialFunctions.hpp>
#include <limits>
#include <stdexcept>

namespace dynsyslab::core::special {
namespace {
constexpr double pi = 3.1415926535897932384626433832795;
void finite(double x) {
    if (!std::isfinite(x))
        throw std::invalid_argument("Special-function arguments must be finite");
}
} // namespace
JacobiValues jacobi(double u, double m) {
    finite(u);
    finite(m);
    if (m < 0) {
        const double a = std::sqrt(1 - m);
        const auto j = jacobi(a * u, m / (m - 1));
        return {j.sn / (a * j.dn), j.cn / j.dn, 1 / j.dn};
    }
    if (m > 1) {
        const double a = std::sqrt(m);
        const auto j = jacobi(a * u, 1 / m);
        return {j.sn / a, j.dn, j.cn};
    }
    JacobiValues j{};
    j.sn = boost::math::jacobi_elliptic(std::sqrt(m), u, &j.cn, &j.dn);
    return j;
}
double ellipticK(double m) {
    finite(m);
    if (m < 0 || m >= 1)
        throw std::domain_error("ellipticK requires 0 <= m < 1");
    return boost::math::ellint_1(std::sqrt(m));
}
double ellipticF(double phi, double m) {
    finite(phi);
    finite(m);
    if (m < 0 || m >= 1)
        throw std::domain_error("ellipticF requires 0 <= m < 1");
    return boost::math::ellint_1(std::sqrt(m), phi);
}
double jacobiAmplitude(double u, double m) {
    if (m < 0 || m >= 1)
        throw std::domain_error("Unwrapped Jacobi amplitude requires 0 <= m < 1");
    const double k = ellipticK(m), turns = std::floor((u + k) / (2 * k)), r = u - 2 * k * turns;
    const auto j = jacobi(r, m);
    return std::atan2(j.sn, j.cn) + turns * pi;
}
double theta(unsigned index, double z, double q) {
    finite(z);
    finite(q);
    if (q < 0 || q >= 1)
        throw std::domain_error("Theta nome must satisfy 0 <= q < 1");
    if (q == 0) {
        if (index < 1 || index > 4)
            throw std::invalid_argument("Theta index must be 1..4");
        return index < 3 ? 0 : 1;
    }
    switch (index) {
    case 1:
        return boost::math::jacobi_theta1(z, q);
    case 2:
        return boost::math::jacobi_theta2(z, q);
    case 3:
        return boost::math::jacobi_theta3(z, q);
    case 4:
        return boost::math::jacobi_theta4(z, q);
    default:
        throw std::invalid_argument("Theta index must be 1..4");
    }
}
WeierstrassValues weierstrass(double z, double g2, double g3) {
    finite(z);
    finite(g2);
    finite(g3);
    const double inf = std::numeric_limits<double>::infinity();
    if (g2 == 0 && g3 == 0) {
        if (z == 0)
            throw std::domain_error("Weierstrass pole");
        return {1 / (z * z), -2 / (z * z * z), inf};
    }
    // Normalize invariants before solving the cubic to avoid overflow/underflow.
    const double scale = std::max(std::sqrt(std::abs(g2)), std::cbrt(std::abs(g3)));
    const double a = g2 / (scale * scale), b = (g3 / scale) / scale / scale;
    const double discriminant = a * a * a - 27 * b * b;
    double w = 0, wp = 0, period = inf;
    double u = z * std::sqrt(scale);
    if (std::abs(discriminant) <= 32 * std::numeric_limits<double>::epsilon() *
                                      std::max(std::abs(a * a * a), 27 * b * b) &&
        a >= 0) {
        const double e = std::sqrt(a / 12), root = std::sqrt(3 * e);
        if (b > 0) {
            period = pi / root;
            u = std::remainder(u, period);
            const double s = std::sin(root * u);
            if (s == 0)
                throw std::domain_error("Weierstrass pole");
            w = -e + 3 * e / (s * s);
            wp = -6 * e * root * std::cos(root * u) / (s * s * s);
        } else {
            const double s = std::sinh(root * u);
            if (s == 0)
                throw std::domain_error("Weierstrass pole");
            w = e + 3 * e / (s * s);
            wp = -6 * e * root * std::cosh(root * u) / (s * s * s);
            if (!std::isfinite(s)) {
                w = e;
                wp = 0;
            }
        }
    } else if (discriminant > 0) {
        const double radius = std::sqrt(a / 3),
                     angle = std::acos(std::clamp(3 * b / a * std::sqrt(3 / a), -1., 1.)) / 3;
        std::array<double, 3> roots{radius * std::cos(angle), radius * std::cos(angle - 2 * pi / 3),
                                    radius * std::cos(angle - 4 * pi / 3)};
        std::sort(roots.begin(), roots.end());
        const double d = roots[2] - roots[0], m = (roots[1] - roots[0]) / d, root = std::sqrt(d);
        period = 2 * ellipticK(m) / root;
        u = std::remainder(u, period);
        const auto j = jacobi(root * u, m);
        if (j.sn == 0)
            throw std::domain_error("Weierstrass pole");
        w = roots[0] + d / (j.sn * j.sn);
        wp = -2 * d * root * j.cn * j.dn / (j.sn * j.sn * j.sn);
    } else {
        const double d = std::sqrt(std::max(0., b * b / 64 - a * a * a / 1728));
        const double e = std::cbrt(b / 8 + d) + std::cbrt(b / 8 - d),
                     h = std::sqrt(3 * e * e - a / 4), root = std::sqrt(h);
        const double m =
            std::clamp(.5 - 3 * e / (4 * h), 0., 1. - std::numeric_limits<double>::epsilon());
        period = 2 * ellipticK(m) / root;
        u = std::remainder(u, period);
        const auto j = jacobi(2 * root * u, m);
        double den = 1 - j.cn;
        if (j.cn > 0)
            den = j.sn * j.sn / (1 + j.cn);
        if (den == 0)
            throw std::domain_error("Weierstrass pole");
        w = e + h * (1 + j.cn) / den;
        wp = -4 * h * root * j.sn * j.dn / (den * den);
    }
    return {scale * w, scale * std::sqrt(scale) * wp, period / std::sqrt(scale)};
}
} // namespace dynsyslab::core::special
