#include "AnalyticSystems.hpp"
#include "HodgkinHuxleyRates.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
namespace dynsyslab::core::detail {
std::vector<AnalyticSolutionFamily> neuralLimitAnalyticSolutions(const SystemConfiguration &c) {
    if (c.identifier != "hodgkin_huxley")
        return {};
    AnalyticSolutionDescriptor d{
        "hh-constant-voltage-gates",
        "Analytic constant-voltage gating",
        R"(V=V_0,\quad q=q_\infty+(q_0-q_\infty)e^{-1000(\alpha_q+\beta_q)(t-t_0)})",
        "Positive capacitance; zero driver; voltage-dependent channel current coefficients vanish; "
        "bias equals leak current at V0",
        "Exact exponentials on a native invariant set",
        {"https://doi.org/10.1113/jphysiol.1952.sp004764"}};
    return {{d, [c](const auto &x, double t0) -> AnalyticTrajectory {
                 const auto &p = c.parameters;
                 if (c.driver.amplitude != 0 || p.at("capacitance") <= 0 ||
                     (p.at("g_na") != 0 && x[0] != p.at("e_na")) ||
                     (p.at("g_k") != 0 && x[0] != p.at("e_k")) ||
                     p.at("bias_current") != p.at("g_l") * (x[0] - p.at("e_l")))
                     return {{},
                             "Native free-voltage current balance does not preserve constant "
                             "voltage for arbitrary gates"};
                 const auto r = hodgkinHuxleyRates(x[0]);
                 const AnalyticState rate{r.am + r.bm, r.ah + r.bh, r.an + r.bn},
                     steady{r.am / rate[0], r.ah / rate[1], r.an / rate[2]};
                 if (!std::all_of(rate.begin(), rate.end(),
                                  [](double v) { return std::isfinite(v) && v > 0; }))
                     return {{}, "Gating rates are not finite and positive at this voltage"};
                 return {[=](double t) {
                             if (t == t0)
                                 return x;
                             auto y = x;
                             for (std::size_t j = 0; j < 3; ++j)
                                 y[j + 1] = steady[j] + (x[j + 1] - steady[j]) *
                                                            std::exp(-1000 * rate[j] * (t - t0));
                             return y;
                         },
                         {}};
             }}};
}
std::vector<AnalyticSolutionFamily> springLimitAnalyticSolutions(const SystemConfiguration &c) {
    if (c.identifier != "spring_coupled_pendulums")
        return {};
    AnalyticSolutionDescriptor d{
        "spring-pendulum-synchronous",
        "Analytic synchronous pendulums",
        R"(\theta_1=\theta_2=\theta,\quad \theta''+(g/l)\sin\theta=0)",
        "Equal positive lengths; spring rest length equals positive pivot separation; zero "
        "damping/input; matching initial angles/velocities",
        "Full nonlinear Jacobi pendulum",
        {"https://dlmf.nist.gov/22.19"}};
    return {{d, [c](const auto &x, double t0) -> AnalyticTrajectory {
                 const auto &p = c.parameters;
                 if (c.driver.amplitude != 0 || p.at("mass1_kg") <= 0 || p.at("mass2_kg") <= 0 ||
                     p.at("length1_m") <= 0 || p.at("length1_m") != p.at("length2_m") ||
                     p.at("pivot_separation_m") <= 1e-12 ||
                     p.at("spring_rest_length_m") != p.at("pivot_separation_m") ||
                     p.at("damping1_n_m_s") != 0 || p.at("damping2_n_m_s") != 0 || x[0] != x[2] ||
                     x[1] != x[3])
                     return {{},
                             "Requires the native conservative synchronous unstretched-spring "
                             "invariant set"};
                 auto orbit =
                     pendulumTrajectory(p.at("gravity_m_s2") / p.at("length1_m"), {x[0], x[1]}, t0);
                 if (!orbit)
                     return orbit;
                 return {[=](double t) {
                             auto a = orbit.evaluate(t);
                             return AnalyticState{a[0], a[1], a[0], a[1]};
                         },
                         {}};
             }}};
}
std::vector<AnalyticSolutionFamily> bubbleLimitAnalyticSolutions(const SystemConfiguration &c) {
    const bool shell = c.identifier == "marmottant";
    if (!shell && c.identifier != "reddy_szeri")
        return {};
    AnalyticSolutionDescriptor d;
    d.identifier = shell ? "marmottant-stationary" : "reddy-szeri-stationary-radius-shapes";
    d.label =
        shell ? "Analytic stationary coated bubble" : "Analytic stationary-radius shape modes";
    d.formula = shell ? R"(r=r_*,\ v=0)"
                      : R"(r=r_*,\ z=z_0,\quad a_n''+\tfrac{Eu C(n-1)(n+1)(n+2)}{2r_*^3}a_n=0)";
    d.conditions = "Zero driver; finite state; verified native radial pressure balance; zero "
                   "radial/translation velocity and effective acoustic forcing";
    d.evaluation = "Exact invariant-set reduction; pressure balance verified at roundoff scale";
    d.references = {shell ? "https://doi.org/10.1121/1.2109427"
                          : "https://doi.org/10.1063/1.1483840"};
    return {
        {d, [c, shell](const auto &x, double t0) -> AnalyticTrajectory {
             if (x.size() != (shell ? 2 : 12) || !std::isfinite(t0) ||
                 !std::all_of(x.begin(), x.end(), [](double a) { return std::isfinite(a); }) ||
                 x[0] <= 0 || x[1] != 0 || c.driver.amplitude != 0)
                 return {{}, "Requires a finite stationary-radius initial state and zero driver"};
             if (!shell &&
                 (x[3] != 0 || (c.parameters.at("translation_enabled") != 0 &&
                                c.dimensionlessGroups.at("angular_frequency") != 0 &&
                                (c.dimensionlessGroups.at("pressure_amplitude") != 0 ||
                                 c.dimensionlessGroups.at("pressure_derivative_amplitude") != 0))))
                 return {{}, "Requires zero translation and zero effective acoustic forcing"};
             // Fresh reset preserves native shell history; no numerical integration is performed.
             auto native = createSystem(c);
             std::vector<double> field(x.size());
             if (native->invalid(t0, x.data()))
                 return {{}, "Initial state is outside the native domain"};
             native->rhs(t0, x.data(), field.data());
             const double eu = c.dimensionlessGroups.at("euler");
             const double scale = std::max(1., std::abs(eu / x[0]));
             if (!std::isfinite(field[1]) ||
                 std::abs(field[1]) > 64 * std::numeric_limits<double>::epsilon() * scale)
                 return {{},
                         "Initial radius does not satisfy the native pressure balance to roundoff"};
             if (shell)
                 return {[x](double) { return x; }, {}};
             std::vector<std::function<AnalyticState(double)>> modes;
             for (int n = 2; n <= 5; ++n) {
                 const double stiffness = .5 * eu * c.dimensionlessGroups.at("capillarity") *
                                          (n - 1) * (n + 1) * (n + 2) / std::pow(x[0], 3);
                 modes.push_back(linearOscillator(0, stiffness, {x[2 * n], x[2 * n + 1]}, t0));
             }
             return {[=](double t) {
                         auto y = x;
                         for (std::size_t j = 0; j < modes.size(); ++j) {
                             const auto a = modes[j](t);
                             y[4 + 2 * j] = a[0];
                             y[5 + 2 * j] = a[1];
                         }
                         return y;
                     },
                     {}};
         }}};
}
} // namespace dynsyslab::core::detail
