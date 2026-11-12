#include "AnalyticSystems.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace dynsyslab::core::detail {
namespace {
constexpr long double pi = 3.141592653589793238462643383279502884L;
// Projective linearization of the Adler equation, preserving the lifted phase.
// (u,v)'=[[K/2,-s],[s,-K/2]](u,v), delta=2 arg(u+iv).
double phaseIncrement(double delta0, double spread, double coupling, double duration) {
    const long double half = std::remainder(static_cast<long double>(delta0) / 2, pi),
                      u = std::cos(half), v = std::sin(half), s = spread,
                      k = static_cast<long double>(coupling) / 2, au = k * u - s * v,
                      av = s * u - k * v;
    const long double discriminant = (std::abs(k) - std::abs(s)) * (std::abs(k) + std::abs(s));
    long double a = 0, b = 0, turns = 0, t = duration;
    if (discriminant < 0) {
        const long double frequency = std::sqrt(-discriminant), period = pi / frequency;
        turns = std::round(t / period);
        if (std::abs(turns) > 1e12L)
            throw std::runtime_error("Adler phase winding exceeds reliable range");
        t -= turns * period;
        const long double co = std::cos(frequency * t), si = std::sin(frequency * t) / frequency;
        a = co * u + si * au;
        b = co * v + si * av;
        turns *= s > 0 ? 1 : -1;
    } else if (discriminant == 0) {
        const long double scale = std::max(1.L, std::abs(t));
        a = u / scale + (t / scale) * au;
        b = v / scale + (t / scale) * av;
    } else {
        const long double rate = std::sqrt(discriminant), z = rate * t;
        if (std::abs(z) < 20) {
            a = std::cosh(z) * u + std::sinh(z) / rate * au;
            b = std::cosh(z) * v + std::sinh(z) / rate * av;
        } else {
            // Scale each eigenspace by its actual coefficient, including exact zero.
            // This retains a pure decaying eigenvector without cosh-sinh cancellation.
            const long double pu = (u + au / rate) / 2, pv = (v + av / rate) / 2,
                              mu = (u - au / rate) / 2, mv = (v - av / rate) / 2,
                              np = std::hypot(pu, pv), nm = std::hypot(mu, mv);
            const long double lp = np == 0 ? -std::numeric_limits<long double>::infinity()
                                           : z + std::log(np),
                              lm = nm == 0 ? -std::numeric_limits<long double>::infinity()
                                           : -z + std::log(nm),
                              top = std::max(lp, lm);
            if (!std::isfinite(top))
                throw std::runtime_error("Adler exponential exceeds reliable range");
            const long double wp = np == 0 ? 0 : std::exp(lp - top) / np,
                              wm = nm == 0 ? 0 : std::exp(lm - top) / nm;
            a = wp * pu + wm * mu;
            b = wp * pv + wm * mv;
        }
    }
    if (!std::isfinite(a) || !std::isfinite(b) || (a == 0 && b == 0))
        throw std::runtime_error("Adler projective coordinates are not resolvable");
    return static_cast<double>(2 * (std::atan2(u * b - v * a, u * a + v * b) + pi * turns));
}
} // namespace
std::vector<AnalyticSolutionFamily> kuramotoPairAnalyticSolutions(const SystemConfiguration &c) {
    if (c.identifier != "kuramoto")
        return {};
    AnalyticSolutionDescriptor d{
        "kuramoto-two-oscillators",
        "Analytic two-oscillator phases",
        R"(\delta'=2s-K\sin\delta,\quad \theta_1+\theta_2=\mathrm{constant})",
        "Exactly two oscillators; zero external driver; finite coupling/frequency spread",
        "Projective matrix exponential with exact phase-slip winding",
        {"https://doi.org/10.4171/EMSS/17"}};
    return {{d, [c](const auto &x, double t0) -> AnalyticTrajectory {
                 const double coupling = c.parameters.at("coupling"),
                              spread = c.parameters.at("frequency_spread");
                 if (x.size() != 2 || c.driver.amplitude != 0 || !std::isfinite(coupling) ||
                     !std::isfinite(spread))
                     return {{}, "Requires the unforced native two-oscillator array"};
                 return {[=](double t) {
                             if (t == t0)
                                 return x;
                             const double increment =
                                 phaseIncrement(x[1] - x[0], spread, coupling, t - t0) / 2;
                             return AnalyticState{x[0] - increment, x[1] + increment};
                         },
                         {}};
             }}};
}
} // namespace dynsyslab::core::detail
