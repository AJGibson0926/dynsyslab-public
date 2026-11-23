#include "AnalyticSystems.hpp"
#include <algorithm>
#include <cmath>
#include <complex>
#include <dynsyslab/core/SpecialFunctions.hpp>
#include <limits>
#include <memory>
#include <stdexcept>

namespace dynsyslab::core::detail {
namespace {
constexpr double pi = 3.1415926535897932384626433832795;
double p(const SystemConfiguration &c, const char *key) { return c.parameters.at(key); }
bool free(const SystemConfiguration &c) { return c.driver.amplitude == 0; }
AnalyticTrajectory no(const char *reason) { return {{}, reason}; }
AnalyticSolutionDescriptor descriptor(const char *id, const char *label, const char *formula,
                                      const char *conditions, const char *reference) {
    return {id, label, formula, conditions, "closed-form special-function evaluation", {reference}};
}
} // namespace
std::function<AnalyticState(double)> linearOscillator(double damping, double stiffness,
                                                      const AnalyticState &initial, double t0,
                                                      SinusoidalDriver driver, double gain) {
    const std::complex<double> root =
        std::sqrt(std::complex<double>(damping * damping / 4 - stiffness));
    auto homogeneous = [=](double t, double x, double v) {
        const auto a = -damping / 2.;
        const auto z = root * t;
        const auto s =
            std::abs(z) < 1e-6 ? t * (1. + z * z / 6. + z * z * z * z / 120.) : std::sinh(z) / root;
        const auto c = std::cosh(z);
        const double decay = std::exp(a * t);
        return AnalyticState{(decay * (x * c + (v - a * x) * s)).real(),
                             (decay * (v * c + (a * v - stiffness * x) * s)).real()};
    };
    const double omega = 2 * pi * driver.frequencyHz;
    const std::complex<double> den(stiffness - omega * omega, damping * omega), i(0, 1);
    if (driver.amplitude == 0)
        return [=](double t) { return homogeneous(t - t0, initial[0], initial[1]); };
    const double amplitude = gain * driver.amplitude;
    if (std::abs(den) > 1e-14 * (1 + std::abs(stiffness) + omega * omega)) {
        const auto coefficient = amplitude / den;
        auto particular = [=](double t) {
            const auto value = coefficient * std::exp(i * (omega * t + driver.phaseRadians));
            return AnalyticState{value.real(), (i * omega * value).real()};
        };
        const auto y0 = particular(t0);
        return [=](double t) {
            auto y = homogeneous(t - t0, initial[0] - y0[0], initial[1] - y0[1]);
            const auto q = particular(t);
            y[0] += q[0];
            y[1] += q[1];
            return y;
        };
    }
    // Exact undamped resonance, including zero-frequency/free-particle forcing.
    auto particular = [=](double t) {
        const double phase = omega * t + driver.phaseRadians;
        if (omega == 0)
            return AnalyticState{.5 * amplitude * std::cos(phase) * t * t,
                                 amplitude * std::cos(phase) * t};
        return AnalyticState{amplitude * t * std::sin(phase) / (2 * omega),
                             amplitude * (std::sin(phase) + omega * t * std::cos(phase)) /
                                 (2 * omega)};
    };
    const auto y0 = particular(t0);
    return [=](double t) {
        auto y = homogeneous(t - t0, initial[0] - y0[0], initial[1] - y0[1]);
        const auto q = particular(t);
        y[0] += q[0];
        y[1] += q[1];
        return y;
    };
}
AnalyticTrajectory pendulumTrajectory(double w2, const AnalyticState &x, double t0) {
    if (x.size() != 2)
        return no("Expected angle and angular velocity");
    if (w2 == 0)
        return {linearOscillator(0, 0, x, t0), {}};
    const double shift = w2 < 0 ? pi : 0, w = std::sqrt(std::abs(w2));
    const double turns = std::floor((x[0] - shift + pi) / (2 * pi)),
                 offset = shift + 2 * pi * turns, a = x[0] - offset;
    const double m = (.5 * x[1] * x[1] + w * w * (1 - std::cos(a))) / (2 * w * w);
    if (x[1] == 0 && (a == 0 || std::abs(std::abs(a) - pi) < 1e-15))
        return {[=](double) { return x; }, {}};
    if (m < 1) {
        if (m == 0)
            return {[=](double) { return x; }, {}};
        const double k = std::sqrt(m), sn = std::clamp(std::sin(a / 2) / k, -1., 1.),
                     cn = x[1] / (2 * k * w);
        const double u0 = special::ellipticF(std::atan2(sn, cn), m);
        return {[=](double t) {
                    const auto j = special::jacobi(u0 + w * (t - t0), m);
                    return AnalyticState{offset + 2 * std::asin(std::clamp(k * j.sn, -1., 1.)),
                                         2 * k * w * j.cn};
                },
                {},4*special::ellipticK(m)/w};
    }
    const double sign = x[1] < 0 ? -1. : 1.;
    if (std::abs(m - 1) <= 4 * std::numeric_limits<double>::epsilon()) {
        const double u0 = std::atanh(std::sin(a / 2));
        return {[=](double t) {
                    const double u = u0 + sign * w * (t - t0);
                    return AnalyticState{offset + 2 * std::asin(std::tanh(u)),
                                         sign * 2 * w / std::cosh(u)};
                },
                {}};
    }
    const double k = 1 / m, frequency = w * std::sqrt(m), u0 = special::ellipticF(a / 2, k);
    return {[=](double t) {
                const double u = u0 + sign * frequency * (t - t0);
                const auto j = special::jacobi(u, k);
                return AnalyticState{offset + 2 * special::jacobiAmplitude(u, k),
                                     sign * 2 * frequency * j.dn};
            },
            {}};
}
std::vector<AnalyticSolutionFamily> oscillatorAnalyticSolutions(const SystemConfiguration &c) {
    std::vector<AnalyticSolutionFamily> result;
    if (c.identifier == "harmonic_oscillator") {
        result.push_back(
            {descriptor(
                 "linear-oscillator", "Analytic oscillator",
                 R"(x''+2\gamma x'+\omega_0^2x=F\cos(\Omega t+\phi))",
                 "Positive mass; constant coefficients; free motion or the system's sinusoidal "
                 "driver",
                 "https://ocw.mit.edu/courses/8-03sc-physics-iii-vibrations-and-waves-fall-2016/"
                 "pages/part-i-mechanical-vibrations-and-waves/"),
             [c](const auto &x, double t0) {
                 if (x.size() != 2)
                     return no("Expected position and velocity");
                 return AnalyticTrajectory{linearOscillator(p(c, "damping_kg_s") / p(c, "mass_kg"),
                                                            p(c, "stiffness_n_m") / p(c, "mass_kg"),
                                                            x, t0, c.driver, 1 / p(c, "mass_kg")),
                                           {},p(c,"damping_kg_s")==0&&free(c)&&p(c,"stiffness_n_m")>0
                                              &&(x[0]!=0||x[1]!=0)?2*pi*std::sqrt(p(c,"mass_kg")/p(c,"stiffness_n_m")):0};
             }});
    } else if (c.identifier == "duffing") {
        result.push_back({descriptor("duffing-elliptic", "Analytic Duffing",
                                     R"(x'^2=2E+\alpha x^2-\tfrac{\beta}{2}x^4)",
                                     "delta=0; no forcing; arbitrary finite ICs in bounded wells, "
                                     "including the separatrix",
                                     "https://dlmf.nist.gov/22.19"),
                          [c](const auto &x, double t0) {
                              if (x.size() != 2)
                                  return no("Expected position and velocity");
                              if (p(c, "delta") != 0 || !free(c))
                                  return no("Undamped, unforced Duffing required");
                              const double a = p(c, "alpha"), b = p(c, "beta");
                              if (b == 0)
                                  return AnalyticTrajectory{linearOscillator(0, -a, x, t0), {}};
                              const double e = .5 * x[1] * x[1] - .5 * a * x[0] * x[0] +
                                               .25 * b * std::pow(x[0], 4);
                              if (b < 0 && !(a < 0 && std::abs(x[0]) <= std::sqrt(a / b) &&
                                             e <= -a * a / (4 * b)))
                                  return no("Escaping softening branch is not covered; no "
                                            "continuation through a finite-time pole");
                              auto orbit = std::make_shared<special::QuarticOrbit>(
                                  std::array<double, 5>{2 * e, 0, a, 0, -b / 2}, x[0], x[1]);
                              return AnalyticTrajectory{[orbit, t0](double t) {
                                                            const auto q = orbit->evaluate(t - t0);
                                                            return AnalyticState{q[0], q[1]};
                                                        },
                                                        {},std::isfinite(orbit->realPeriod())?orbit->realPeriod():0};
                          }});
        result.push_back({descriptor("duffing-linear", "Analytic linear Duffing limit",
                                     R"(x''+\delta x'-\alpha x=g u(t))",
                                     "beta=0; free motion or sinusoidal system driver",
                                     "https://dlmf.nist.gov/1.13"),
                          [c](const auto &x, double t0) {
                              if (p(c, "beta") != 0)
                                  return no("beta must be zero");
                              return AnalyticTrajectory{
                                  linearOscillator(p(c, "delta"), -p(c, "alpha"), x, t0, c.driver,
                                                   p(c, "input_gain")),
                                  {}};
                          }});
    } else if (c.identifier == "nonlinear_pendulum") {
        result.push_back({descriptor("pendulum-jacobi", "Analytic pendulum",
                                     R"(\theta=2\arcsin(k\operatorname{sn}(u|k^2)))",
                                     "No damping or forcing; librations, rotations and separatrix; "
                                     "unwrapped angle preserved",
                                     "https://dlmf.nist.gov/22.19#i"),
                          [c](const auto &x, double t0) {
                              if (p(c, "damping_n_m_s") != 0 || !free(c))
                                  return no("Zero damping and forcing required");
                              return pendulumTrajectory(p(c, "gravity_m_s2") / p(c, "length_m"), x,
                                                        t0);
                          }});
        result.push_back(
            {descriptor("pendulum-free-rotor", "Analytic free rotor", R"(\theta''+\gamma\theta'=0)",
                        "gravity=0; no forcing", "https://dlmf.nist.gov/1.13"),
             [c](const auto &x, double t0) {
                 if (p(c, "gravity_m_s2") != 0 || !free(c))
                     return no("Zero gravity and forcing required");
                 return AnalyticTrajectory{
                     linearOscillator(p(c, "damping_n_m_s") /
                                          (p(c, "mass_kg") * std::pow(p(c, "length_m"), 2)),
                                      0, x, t0),
                     {}};
             }});
    } else if (c.identifier == "van_der_pol") {
        result.push_back({descriptor("vdp-harmonic-limit", "Analytic harmonic limit",
                                     R"(x(t)=x_0\cos t+v_0\sin t)",
                                     "mu=0; free motion or sinusoidal system driver",
                                     "https://www.scholarpedia.org/article/Van_der_Pol_oscillator"),
                          [c](const auto &x, double t0) {
                              if (p(c, "mu") != 0)
                                  return no("mu must be zero");
                              return AnalyticTrajectory{linearOscillator(0, 1, x, t0, c.driver),
                                                        {}};
                          }});
    }
    return result;
}
} // namespace dynsyslab::core::detail
