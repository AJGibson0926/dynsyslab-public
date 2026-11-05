#include "AnalyticSystems.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <dynsyslab/core/SpecialFunctions.hpp>
#include <limits>
namespace dynsyslab::core::detail {
namespace {
double potential(double q) { return q * q * (.5 - q / 3); }
double turning(double energy, double lo, double hi, bool increasing) {
    for (int k = 0; k < 80; ++k) {
        const double mid = (lo + hi) / 2;
        if ((potential(mid) < energy) == increasing)
            lo = mid;
        else
            hi = mid;
    }
    return (lo + hi) / 2;
}
AnalyticTrajectory straightLine(const SystemConfiguration &c, const AnalyticState &x, double t0) {
    const double lam = c.parameters.at("lam");
    if (c.driver.amplitude != 0 || lam == 0)
        return {{}, "Requires nonzero lam and no forcing"};
    const std::array<std::array<double, 2>, 3> axes{
        {{0, 1}, {std::sqrt(3.) / 2, -.5}, {-std::sqrt(3.) / 2, -.5}}};
    double ex = 0, ey = 0;
    bool found = false;
    for (const auto &e : axes) {
        const double q = std::hypot(x[0], x[1]), v = std::hypot(x[2], x[3]),
                     tol = 32 * std::numeric_limits<double>::epsilon();
        if (std::abs(e[1] * x[0] - e[0] * x[1]) <= tol * q &&
            std::abs(e[1] * x[2] - e[0] * x[3]) <= tol * v) {
            ex = e[0];
            ey = e[1];
            found = true;
            break;
        }
    }
    if (!found)
        return {{}, "Initial position and momentum must lie on the same exact symmetry axis"};
    const double q = lam * (ex * x[0] + ey * x[1]), v = lam * (ex * x[2] + ey * x[3]);
    auto lift = [=](double y, double w) {
        return AnalyticState{ex * y / lam, ey * y / lam, ex * w / lam, ey * w / lam};
    };
    if (v == 0 && (q == 0 || q == 1))
        return {[x](double) { return x; }, {}};
    const double energy = .5 * v * v + potential(q), barrier = 1. / 6;
    if (!std::isfinite(energy) || q >= 1 || energy <= 0 || energy > barrier)
        return {
            {},
            "Requires the bounded well or its homoclinic boundary; escaping branches are excluded"};
    if (energy == barrier) {
        const double ratio = v / (1 - q);
        if (std::abs(ratio) >= 1)
            return {{}, "Degenerate separatrix phase is not representable"};
        const double z0 = std::atanh(ratio);
        return {[=](double t) {
                    const double z = z0 + (t - t0) / 2, a = std::exp(-2 * std::abs(z)),
                                 sech2 = 4 * a / ((1 + a) * (1 + a));
                    return lift(1 - 1.5 * sech2, 1.5 * sech2 * std::tanh(z));
                },
                {}};
    }
    const double q1 = turning(energy, -.5, 0, false), q2 = turning(energy, 0, 1, true),
                 q3 = turning(energy, 1, 1.5, false);
    const double span = q2 - q1, m = span / (q3 - q1), rate = std::sqrt((q3 - q1) / 6),
                 K = special::ellipticK(m);
    double phase = special::ellipticF(std::asin(std::sqrt(std::clamp((q - q1) / span, 0., 1.))), m);
    if (v < 0)
        phase = 2 * K - phase;
    return {[=](double t) {
                if (t == t0)
                    return x;
                const auto j = special::jacobi(std::remainder(rate * (t - t0) + phase, 2 * K), m);
                return lift(q1 + span * j.sn * j.sn, 2 * rate * span * j.sn * j.cn * j.dn);
            },
            {}};
}
} // namespace
std::vector<AnalyticSolutionFamily> henonHeilesAnalyticSolutions(const SystemConfiguration &c) {
    return {{{"henon-heiles-straight-line",
              "Analytic symmetry-axis orbit",
              R"(q=q_1+(q_2-q_1)\operatorname{sn}^2(a(t-t_0)+z_0\mid m),\quad q=\lambda s)",
              "No forcing; position and momentum on one of the three symmetry axes; bounded well "
              "0<=lambda^2 H<=1/6 or a stationary saddle",
              "Jacobi elliptic evaluation; elementary homoclinic boundary",
              {"https://epub.uni-regensburg.de/1660/1/brack11.pdf"}},
             [c](const auto &x, double t0) { return straightLine(c, x, t0); }}};
}
} // namespace dynsyslab::core::detail
