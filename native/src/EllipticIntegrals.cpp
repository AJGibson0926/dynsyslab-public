// Carlson duplication for R_F, R_D, R_C and R_J (Carlson 1995, Numer. Algorithms 10:13-26; DLMF
// 19.36.1, 19.36.2) and the Legendre integrals in the parameter convention (DLMF 19.25.5, 19.25.7,
// 19.25.9, 19.25.14).
#include <dynsyslab/core/EllipticIntegrals.hpp>

#include <algorithm>
#include <cmath>
#include <limits>

namespace dynsyslab::core {
namespace {
constexpr double nan = std::numeric_limits<double>::quiet_NaN();
}

double carlsonRF(double x, double y, double z) {
    if (!(x >= 0 && y >= 0 && z >= 0) || (x + y == 0) || (x + z == 0) || (y + z == 0)) return nan;
    for (int iteration = 0; iteration < 200; ++iteration) {
        const double mu = (x + y + z) / 3;
        const double dx = 1 - x / mu, dy = 1 - y / mu, dz = 1 - z / mu;
        if (std::max({std::abs(dx), std::abs(dy), std::abs(dz)}) < 1e-4) {
            const double e2 = dx * dy - dz * dz, e3 = dx * dy * dz;
            return (1 - e2 / 10 + e3 / 14 + e2 * e2 / 24 - 3 * e2 * e3 / 44) / std::sqrt(mu);
        }
        const double sx = std::sqrt(x), sy = std::sqrt(y), sz = std::sqrt(z);
        const double lambda = sx * (sy + sz) + sy * sz;
        x = (x + lambda) / 4; y = (y + lambda) / 4; z = (z + lambda) / 4;
    }
    return nan;
}

double carlsonRD(double x, double y, double z) {
    if (!(x >= 0 && y >= 0 && z > 0) || x + y == 0) return nan;
    double sum = 0, factor = 1;
    for (int iteration = 0; iteration < 200; ++iteration) {
        const double mu = (x + y + 3 * z) / 5;
        const double dx = 1 - x / mu, dy = 1 - y / mu, dz = 1 - z / mu;
        if (std::max({std::abs(dx), std::abs(dy), std::abs(dz)}) < 1e-4) {
            const double ea = dx * dy, eb = dz * dz, ec = ea - eb, ed = ea - 6 * eb, ee = ed + ec + ec;
            const double series = 1 + ed * (-3.0 / 14 + 9.0 / 88 * ed - 4.5 / 26 * dz * ee)
                                + dz * (ee / 6 + dz * (-9.0 / 22 * ec + 3.0 / 26 * dz * ea));
            return 3 * sum + factor * series / (mu * std::sqrt(mu));
        }
        const double sx = std::sqrt(x), sy = std::sqrt(y), sz = std::sqrt(z);
        const double lambda = sx * (sy + sz) + sy * sz;
        sum += factor / (sz * (z + lambda));
        factor /= 4;
        x = (x + lambda) / 4; y = (y + lambda) / 4; z = (z + lambda) / 4;
    }
    return nan;
}

double carlsonRC(double x, double y) {
    if (!(x >= 0 && y > 0)) return nan;
    for (int iteration = 0; iteration < 200; ++iteration) {
        const double mu = (x + 2 * y) / 3, s = (y - mu) / mu;
        if (std::abs(s) < 1e-4)
            return (1 + s * s * (3.0 / 10 + s * (1.0 / 7 + s * (3.0 / 8 + s * 9.0 / 22)))) / std::sqrt(mu);
        const double lambda = 2 * std::sqrt(x) * std::sqrt(y) + y;
        x = (x + lambda) / 4; y = (y + lambda) / 4;
    }
    return nan;
}

double carlsonRJ(double x, double y, double z, double p) {
    if (!(x >= 0 && y >= 0 && z >= 0 && p > 0) || x + y == 0 || x + z == 0 || y + z == 0) return nan;
    double sum = 0, factor = 1;
    for (int iteration = 0; iteration < 200; ++iteration) {
        const double mu = (x + y + z + 2 * p) / 5;
        const double dx = 1 - x / mu, dy = 1 - y / mu, dz = 1 - z / mu, dp = 1 - p / mu;
        if (std::max({std::abs(dx), std::abs(dy), std::abs(dz), std::abs(dp)}) < 1e-4) {
            // DLMF 19.36.2, with X + Y + Z + 2P = 0.
            const double e2 = dx * dy + dx * dz + dy * dz - 3 * dp * dp, xyz = dx * dy * dz;
            const double e3 = xyz + 2 * e2 * dp + 4 * dp * dp * dp, e4 = (2 * xyz + e2 * dp + 3 * dp * dp * dp) * dp,
                         e5 = xyz * dp * dp;
            const double series = 1 - 3 * e2 / 14 + e3 / 6 + 9 * e2 * e2 / 88 - 3 * e4 / 22 - 9 * e2 * e3 / 52 + 3 * e5 / 26;
            return 3 * sum + factor * series / (mu * std::sqrt(mu));
        }
        const double sx = std::sqrt(x), sy = std::sqrt(y), sz = std::sqrt(z);
        const double lambda = sx * (sy + sz) + sy * sz;
        const double alpha = std::pow(p * (sx + sy + sz) + sx * sy * sz, 2), beta = p * (p + lambda) * (p + lambda);
        sum += factor * carlsonRC(alpha, beta);
        factor /= 4;
        x = (x + lambda) / 4; y = (y + lambda) / 4; z = (z + lambda) / 4; p = (p + lambda) / 4;
    }
    return nan;
}

double ellipticK(double m) { return m < 1 ? carlsonRF(0, 1 - m, 1) : nan; }

double ellipticE(double m) {
    if (!(m <= 1)) return nan;
    if (m == 1) return 1;
    return carlsonRF(0, 1 - m, 1) - m / 3 * carlsonRD(0, 1 - m, 1);
}

double ellipticF(double phi, double m) {
    const double turns = std::round(phi / std::acos(-1.0)), r = phi - turns * std::acos(-1.0);
    const double s = std::sin(r), c = std::cos(r);
    const double base = s * carlsonRF(c * c, 1 - m * s * s, 1);
    return turns == 0 ? base : base + 2 * turns * ellipticK(m);
}

double ellipticEIncomplete(double phi, double m) {
    const double turns = std::round(phi / std::acos(-1.0)), r = phi - turns * std::acos(-1.0);
    const double s = std::sin(r), c = std::cos(r), delta = 1 - m * s * s;
    const double base = s * carlsonRF(c * c, delta, 1) - m / 3 * s * s * s * carlsonRD(c * c, delta, 1);
    return turns == 0 ? base : base + 2 * turns * ellipticE(m);
}

double ellipticPi(double phi, double n, double m) {
    if (!(n < 1)) return nan;
    const double turns = std::round(phi / std::acos(-1.0)), r = phi - turns * std::acos(-1.0);
    const double s = std::sin(r), c = std::cos(r), delta = 1 - m * s * s;
    const double base = s * carlsonRF(c * c, delta, 1) + n / 3 * s * s * s * carlsonRJ(c * c, delta, 1, 1 - n * s * s);
    if (turns == 0) return base;
    const double complete = carlsonRF(0, 1 - m, 1) + n / 3 * carlsonRJ(0, 1 - m, 1, 1 - n);
    return base + 2 * turns * complete;
}

}  // namespace dynsyslab::core
