#include "AnalyticSystems.hpp"
#include <cmath>
#include <dynsyslab/core/SpecialFunctions.hpp>
#include <memory>

namespace dynsyslab::core::detail {
namespace {
constexpr double pi = 3.14159265358979323846;
double p(const SystemConfiguration &c, const char *key) { return c.parameters.at(key); }
AnalyticTrajectory no(const char *reason) { return {{}, reason}; }
AnalyticSolutionFamily
family(const char *id, const char *label, const char *formula, const char *conditions,
       std::function<AnalyticTrajectory(const AnalyticState &, double)> make) {
    return {{id,
             label,
             formula,
             conditions,
             "closed-form evaluation",
             {"https://www.mat.univie.ac.at/~gerald/ftp/book-ode/ode.pdf"}},
            std::move(make)};
}
AnalyticTrajectory affine(std::vector<double> a, AnalyticState b, const AnalyticState &x,
                          double t0) {
    return {[=](double t) { return special::affineExponential(a, b, x, t - t0); }, {}};
}
} // namespace
std::vector<AnalyticSolutionFamily> otherAnalyticSolutions(const SystemDescriptor &d,
                                                           const SystemConfiguration &c) {
    if (d.identifier == "double_pendulum") return doublePendulumAnalyticSolutions(c);
    std::vector<AnalyticSolutionFamily> result;
    if (d.identifier == "lorenz") {
        result.push_back(family(
            "lorenz-invariant-axis", "Analytic invariant axis", R"(x=y=0,\ z=z_0e^{-\beta(t-t_0)})",
            "x0=y0=0; no forcing", [c](const auto &x, double t0) -> AnalyticTrajectory {
                if (c.driver.amplitude != 0 || x[0] != 0 || x[1] != 0)
                    return no("Requires the unforced invariant z axis");
                return {[=](double t) {
                            return AnalyticState{0, 0, x[2] * std::exp(-p(c, "beta") * (t - t0))};
                        },
                        {}};
            }));
        result.push_back(family("lorenz-zero-sigma", "Analytic affine Lorenz limit",
                                R"(x=x_0,\ (y,z)'=A(x_0)(y,z)+b(x_0))", "sigma=0; no forcing",
                                [c](const auto &x, double t0) -> AnalyticTrajectory {
                                    if (c.driver.amplitude != 0 || p(c, "sigma") != 0)
                                        return no("Requires sigma=0 and no forcing");
                                    auto q = affine({-1, -x[0], x[0], -p(c, "beta")},
                                                    {x[0] * p(c, "rho"), 0}, {x[1], x[2]}, t0);
                                    return {[q, x](double t) {
                                                auto y = q.evaluate(t);
                                                return AnalyticState{x[0], y[0], y[1]};
                                            },
                                            {}};
                                }));
    } else if (d.identifier == "rossler") {
        result.push_back(
            family("rossler-invariant-plane", "Analytic invariant plane",
                   R"(z=0,\ (x,y)'=\begin{pmatrix}0&-1\\1&a\end{pmatrix}(x,y))",
                   "b=0, z0=0; no forcing", [c](const auto &x, double t0) -> AnalyticTrajectory {
                       if (c.driver.amplitude != 0 || p(c, "b") != 0 || x[2] != 0)
                           return no("Requires b=0, z0=0 and no forcing");
                       auto q = affine({0, -1, 1, p(c, "a")}, {0, 0}, {x[0], x[1]}, t0);
                       return {[q](double t) {
                                   auto y = q.evaluate(t);
                                   y.push_back(0);
                                   return y;
                               },
                               {}};
                   }));
    } else if (d.identifier == "chua") {
        result.push_back(
            family("chua-linear-limit", "Analytic linear circuit", R"(X(t)=e^{A(t-t_0)}X_0)",
                   "m0=m1; no forcing", [c](const auto &x, double t0) -> AnalyticTrajectory {
                       if (c.driver.amplitude != 0 || p(c, "m0") != p(c, "m1"))
                           return no("Requires equal Chua slopes and no forcing");
                       double a = p(c, "alpha");
                       return affine({-a * (1 + p(c, "m0")), a, 0, 1, -1, 1, 0, -p(c, "beta"), 0},
                                     {0, 0, 0}, x, t0);
                   }));
    } else if (d.identifier == "hindmarsh_rose") {
        result.push_back(family(
            "hr-affine-limit", "Analytic affine neuron limit", R"(X'=AX+b)", "a=b=d=0; no forcing",
            [c](const auto &x, double t0) -> AnalyticTrajectory {
                if (c.driver.amplitude != 0 || p(c, "a") != 0 || p(c, "b") != 0 || p(c, "d") != 0)
                    return no("Requires a=b=d=0 and no forcing");
                double r = p(c, "r"), s = p(c, "s");
                return affine({0, 1, -1, 0, -1, 0, r * s, 0, -r},
                              {p(c, "current"), p(c, "c"), -r * s * p(c, "x_r")}, x, t0);
            }));
    } else if (d.identifier == "henon_heiles") {
        result = henonHeilesAnalyticSolutions(c);
        result.push_back(family("henon-heiles-harmonic", "Analytic harmonic limit",
                                R"(q_j''+q_j=0)", "lam=0; no forcing",
                                [c](const auto &x, double t0) -> AnalyticTrajectory {
                                    if (c.driver.amplitude != 0 || p(c, "lam") != 0)
                                        return no("Requires lam=0 and no forcing");
                                    auto a = linearOscillator(0, 1, {x[0], x[2]}, t0),
                                         b = linearOscillator(0, 1, {x[1], x[3]}, t0);
                                    return {[=](double t) {
                                                auto y = a(t), z = b(t);
                                                return AnalyticState{y[0], z[0], y[1], z[1]};
                                            },
                                            {}};
                                }));
    } else if (d.identifier == "fput") {
        result.push_back(
            family("fput-normal-modes", "Analytic normal modes",
                   R"(q_j=\sum_k Q_k\sin(\pi jk/(N+1)),\ \omega_k=2\sin(\pi k/(2N+2)))",
                   "alpha=beta=0; fixed ends; no forcing",
                   [c](const auto &x, double t0) -> AnalyticTrajectory {
                       if (c.driver.amplitude != 0 || p(c, "alpha") != 0 || p(c, "beta") != 0)
                           return no("Requires the unforced linear FPUT chain");
                       const std::size_t n = x.size() / 2;
                       AnalyticState q(n, 0), v(n, 0), w(n);
                       for (std::size_t k = 0; k < n; ++k) {
                           w[k] = 2 * std::sin(pi * (k + 1) / (2 * (n + 1)));
                           for (std::size_t j = 0; j < n; ++j) {
                               double s = 2 * std::sin(pi * (j + 1) * (k + 1) / (n + 1)) / (n + 1);
                               q[k] += s * x[j];
                               v[k] += s * x[n + j];
                           }
                       }
                       return {[=](double t) {
                                   AnalyticState y(2 * n, 0);
                                   for (std::size_t k = 0; k < n; ++k) {
                                       double co = std::cos(w[k] * (t - t0)),
                                              si = std::sin(w[k] * (t - t0));
                                       for (std::size_t j = 0; j < n; ++j) {
                                           double s = std::sin(pi * (j + 1) * (k + 1) / (n + 1));
                                           y[j] += s * (q[k] * co + v[k] * si / w[k]);
                                           y[n + j] += s * (-w[k] * q[k] * si + v[k] * co);
                                       }
                                   }
                                   return y;
                               },
                               {}};
                   }));
    } else if (d.identifier == "kuramoto") {
        result.push_back(
            family("kuramoto-uncoupled", "Analytic uncoupled phases",
                   R"(\theta_j=\theta_{j0}+\omega_j(t-t_0))", "coupling=0 or one oscillator; no forcing",
                   [c](const auto &x, double t0) -> AnalyticTrajectory {
                       if (c.driver.amplitude != 0 || (x.size() != 1 && p(c, "coupling") != 0))
                           return no("Requires zero coupling and no forcing");
                       return {[=](double t) {
                                   auto y = x;
                                   for (std::size_t j = 0; j < x.size(); ++j) {
                                       double f = x.size() == 1 ? 0 : double(j) / (x.size() - 1);
                                       y[j] += p(c, "frequency_spread") * (2 * f - 1) * (t - t0);
                                   }
                                   return y;
                               },
                               {}};
                   }));
    } else if (d.identifier == "spring_coupled_pendulums") {
        result.push_back(family(
            "decoupled-pendulums", "Analytic decoupled pendulums",
            R"(\theta_j''+(g/l_j)\sin\theta_j=0)", "spring=0; damping1=damping2=0; no forcing",
            [c](const auto &x, double t0) -> AnalyticTrajectory {
                if (c.driver.amplitude != 0 || p(c, "spring_constant_n_m") != 0 ||
                    p(c, "damping1_n_m_s") != 0 || p(c, "damping2_n_m_s") != 0)
                    return no("Requires uncoupled conservative pendulums");
                auto a =
                    pendulumTrajectory(p(c, "gravity_m_s2") / p(c, "length1_m"), {x[0], x[1]}, t0);
                auto b =
                    pendulumTrajectory(p(c, "gravity_m_s2") / p(c, "length2_m"), {x[2], x[3]}, t0);
                return {[=](double t) {
                            auto y = a.evaluate(t), z = b.evaluate(t);
                            y.insert(y.end(), z.begin(), z.end());
                            return y;
                        },
                        {}};
            }));
    }
    auto neural = neuralLimitAnalyticSolutions(c), spring = springLimitAnalyticSolutions(c);
    result.insert(result.end(), neural.begin(), neural.end());
    result.insert(result.end(), spring.begin(), spring.end());
    auto pair = kuramotoPairAnalyticSolutions(c);
    result.insert(result.end(), pair.begin(), pair.end());
    return result;
}
} // namespace dynsyslab::core::detail
