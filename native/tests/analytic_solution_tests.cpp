#include <cmath>
#include <dynsyslab/core/Integrator.hpp>
#include <dynsyslab/core/SpecialFunctions.hpp>
#include <dynsyslab/core/System.hpp>
#include <iostream>
#include <stdexcept>
using namespace dynsyslab::core;
namespace {
void check(bool yes, const std::string &label) {
    if (!yes)
        throw std::runtime_error(label);
}
void near(double a, double b, double tolerance, const std::string &label) {
    if (!std::isfinite(a) || !std::isfinite(b) || std::abs(a - b) > tolerance * (1 + std::abs(b)))
        throw std::runtime_error(label + ": " + std::to_string(a) + " != " + std::to_string(b));
}
void specialFunctions() {
    near(special::ellipticK(.5), 1.8540746773013719184, 2e-15, "elliptic K");
    auto j = special::jacobi(1, .5);
    near(j.sn, .80300182489564388764, 2e-15, "Jacobi sn reference");
    near(j.cn, .59597656767214067402, 2e-15, "Jacobi cn reference");
    near(j.dn, .82316100163159626945, 2e-15, "Jacobi dn reference");
    near(special::theta(3, 0, .1), 1.2002000020000002, 2e-15, "theta series");
    for (double m : {-2., 0., .4, 1., 2.})
        for (double u : {-.3, .1, 1.3, 7.}) {
            auto q = special::jacobi(u, m);
            near(q.sn * q.sn + q.cn * q.cn, 1, 3e-14, "Jacobi circle identity");
            near(q.dn * q.dn + m * q.sn * q.sn, 1, 3e-14, "Jacobi elliptic identity");
        }
    for (auto invariants : std::vector<std::vector<double>>{
             {4, 0}, {0, 1}, {-3, 2}, {12, 8}, {12, -8}, {0, 0}, {1e12, 3e17}}) {
        double g2 = invariants[0], g3 = invariants[1],
               scale = std::max(1., std::pow(std::abs(g2), .25));
        for (double arg : {.15, .3, .7, 1.1}) {
            double z = arg / scale;
            auto w = special::weierstrass(z, g2, g3);
            near(w.derivative * w.derivative, 4 * w.value * w.value * w.value - g2 * w.value - g3,
                 2e-11, "Weierstrass differential identity");
        }
    }
}
void compare(SystemConfiguration c, AnalyticState x, double duration = 12,
             double tolerance = 2e-7) {
    auto system = createSystem(c);
    auto families = system->analyticSolutions();
    check(!families.empty(), c.identifier + " owns families");
    IntegrationOptions options;
    options.initialTime = .371;
    options.finalTime = options.initialTime + duration;
    options.outputStep = .037;
    options.relativeTolerance = 1e-11;
    options.absoluteTolerance = 1e-13;
    options.maximumStep = .012;
    auto numerical = integrate(*system, x, 1, options);
    auto result =
        compareAnalyticSolution(families, x, options.initialTime, numerical.time, numerical.values);
    check(result.status == "complete",
          c.identifier + " complete: " + result.status + " " + result.reason);
    for (std::size_t i = 0; i < result.values.size(); ++i)
        near(result.values[i], numerical.values[i], tolerance,
             c.identifier + " sample " + std::to_string(i));
    auto instance = families.front().instantiate(x, options.initialTime);
    if (instance) {
        const double t = options.initialTime + .271, h = 1e-5;
        auto y = instance.evaluate(t), a = instance.evaluate(t - h), b = instance.evaluate(t + h);
        std::vector<double> rhs(x.size());
        system->rhs(t, y.data(), rhs.data());
        for (std::size_t i = 0; i < x.size(); ++i)
            near((b[i] - a[i]) / (2 * h), rhs[i], 2e-6, c.identifier + " ODE residual");
    }
    auto disabled = compareAnalyticSolution(families, x, options.initialTime, numerical.time,
                                            numerical.values, "off");
    check(disabled.status == "disabled" && disabled.values.empty(), "explicit analytic opt out");
    check(!families.front().instantiate({}, 0), "reject malformed IC");
    std::cout << c.identifier << " " << result.descriptor.identifier << " checked\n";
}
SystemConfiguration config(const char *id, std::map<std::string, double> parameters = {}) {
    SystemConfiguration c;
    c.identifier = id;
    c.parameters = std::move(parameters);
    return c;
}
void oscillators() {
    for (double damping : {0., .3, 2., 4.})
        compare(config("harmonic_oscillator", {{"damping_kg_s", damping}}), {.6, -.4});
    for (double frequency : {.11, 1 / (2 * 3.14159265358979323846)}) {
        auto c = config("harmonic_oscillator");
        c.driver = {.3, frequency, .2};
        compare(c, {.6, -.4});
    }
    for (double alpha : {-1., 0., 1.})
        for (auto x : std::vector<AnalyticState>{{.3, .4}, {1.4, 0}, {-.3, -.4}, {0, 1.5}})
            compare(config("duffing", {{"delta", 0}, {"alpha", alpha}}), x);
    compare(config("duffing", {{"delta", 0}}), {std::sqrt(2.), 0});
    compare(config("duffing", {{"delta", 0}}), {1, 0});
    compare(config("duffing", {{"delta", .2}, {"beta", 0}}), {.3, .2}, 4);
    for (auto x : std::vector<AnalyticState>{{.6, 0}, {.6, -.4}, {7., 8.}, {7., -8.}})
        compare(config("nonlinear_pendulum", {{"damping_n_m_s", 0}, {"gravity_m_s2", 1}}), x, 18);
    compare(config("nonlinear_pendulum", {{"damping_n_m_s", 0}, {"gravity_m_s2", 1}}), {0, 2}, 6);
    compare(config("van_der_pol", {{"mu", 0}}), {.2, .4});
}
void catalogue() {
    compare(config("lorenz"), {0, 0, 8});
    compare(config("lorenz", {{"sigma", 0}}), {.3, 1, 2});
    compare(config("rossler", {{"b", 0}}), {.3, 1, 0});
    compare(config("chua", {{"m0", -.2}, {"m1", -.2}}), {.3, .1, 0}, 3);
    compare(config("hindmarsh_rose", {{"a", 0}, {"b", 0}, {"d", 0}}), {.3, 1, 2});
    compare(config("henon_heiles", {{"lam", 0}}), {.2, .1, .4, 0});
    compare(config("fput", {{"sites", 3}, {"alpha", 0}, {"beta", 0}}), {.2, 0, -.1, .4, .2, 0});
    compare(config("kuramoto", {{"oscillators", 3}, {"coupling", 0}}), {.2, .1, 4});
    compare(config("spring_coupled_pendulums",
                   {{"spring_constant_n_m", 0}, {"damping1_n_m_s", 0}, {"damping2_n_m_s", 0}}),
            {.2, .1, .4, .2});
    for (double kappa : {4. / 3., 1.5}) {
        auto c = config("rpe",
                        {{"kappa", kappa}, {"mu_pa_s", 0}, {"sigma_n_m", kappa == 1.5 ? 0 : .072}});
        c.preset = "air_water_20c";
        compare(c, {1.15, .04}, 6, 2e-7);
    }
    auto damped = createSystem(config("duffing"));
    check(!damped->analyticSolutions().front().instantiate({.2, 0}, 0),
          "damped nonlinear Duffing excluded");
}
} // namespace
int main() {
    try {
        specialFunctions();
        oscillators();
        catalogue();
        std::cout << "analytic solutions: PASS\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
