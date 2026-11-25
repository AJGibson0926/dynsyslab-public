// Carlson duplication for R_F and R_D (Carlson 1995, Numer. Algorithms 10:13-26; DLMF 19.36.1) and
// the Legendre integrals in the parameter convention (DLMF 19.25.5, 19.25.7, 19.25.9).
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

}  // namespace dynsyslab::core
