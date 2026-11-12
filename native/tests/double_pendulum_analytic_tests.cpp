#include <algorithm>
#include <cmath>
#include <dynsyslab/core/Integrator.hpp>
#include <dynsyslab/core/SpecialFunctions.hpp>
#include <dynsyslab/core/System.hpp>
#include <iostream>
#include <stdexcept>
using namespace dynsyslab::core;
namespace {
constexpr double pi = 3.14159265358979323846;
void check(bool ok, const std::string &message) {
    if (!ok)
        throw std::runtime_error(message);
}
void near(double a, double b, double tolerance, const std::string &message) {
    if (!std::isfinite(a) || !std::isfinite(b) || std::abs(a - b) > tolerance * (1 + std::abs(b)))
        throw std::runtime_error(message + ": " + std::to_string(a) + " versus " +
                                 std::to_string(b));
}
SystemConfiguration config(double m1 = 1, double m2 = 1, double l1 = 1, double l2 = 1) {
    SystemConfiguration c;
    c.identifier = "double_pendulum";
    c.parameters = {{"mass1_kg", m1},     {"mass2_kg", m2},    {"length1_m", l1},
                    {"length2_m", l2},    {"gravity_m_s2", 0}, {"damping1_n_m_s", 0},
                    {"damping2_n_m_s", 0}};
    return c;
}
std::pair<double, double> invariants(const SystemConfiguration &c, const AnalyticState &x) {
    const auto &p = c.parameters;
    const double a = (p.at("mass1_kg") + p.at("mass2_kg")) * std::pow(p.at("length1_m"), 2),
                 b = p.at("mass2_kg") * p.at("length1_m") * p.at("length2_m"),
                 cc = p.at("mass2_kg") * std::pow(p.at("length2_m"), 2), co = std::cos(x[0] - x[2]);
    return {.5 * (a * x[1] * x[1] + 2 * b * co * x[1] * x[3] + cc * x[3] * x[3]),
            (a + b * co) * x[1] + (cc + b * co) * x[3]};
}
void compare(const SystemConfiguration &c, const AnalyticState &initial, double duration = 30) {
    std::cout << "Checking DP " << initial[0] << "," << initial[1] << "," << initial[2] << ","
              << initial[3] << std::endl;
    auto system = createSystem(c);
    auto families = system->analyticSolutions();
    check(families.size() == 1, "double pendulum owns its analytic family");
    auto orbit = families[0].instantiate(initial, .371);
    check(bool(orbit), orbit.reason);
    check(orbit.evaluate(.371) == initial, "exact IC identity");
    const auto invariant = invariants(c, initial);
    for (double t : {-5., -.01, .2, 2., 9., 20., 50.}) {
        const auto y = orbit.evaluate(t), a = orbit.evaluate(t - 1e-5),
                   b = orbit.evaluate(t + 1e-5);
        const auto inv = invariants(c, y);
        near(inv.first, invariant.first, 2e-10, "conserved energy");
        near(inv.second, invariant.second, 2e-10, "conserved angular momentum");
        AnalyticState field(4);
        system->rhs(t, y.data(), field.data());
        for (std::size_t j = 0; j < 4; ++j)
            near((b[j] - a[j]) / 2e-5, field[j], 3e-6, "native RHS residual");
    }
    double previous = 1e100;
    for (double tolerance : {1e-7, 1e-10}) {
        IntegrationOptions options;
        options.initialTime = .371;
        options.finalTime = .371 + duration;
        options.outputStep = .083;
        options.relativeTolerance = tolerance;
        options.absoluteTolerance = tolerance * .01;
        options.maximumStep = .08;
        const auto numerical = integrate(*system, initial, 1, options);
        const auto compared =
            compareAnalyticSolution(families, initial, .371, numerical.time, numerical.values);
        check(compared.status == "complete", "complete analytic comparison: " + compared.reason);
        double error = 0;
        for (std::size_t i = 0; i < compared.values.size(); ++i)
            error = std::max(error, std::abs(compared.values[i] - numerical.values[i]));
        std::cout << "  tolerance " << tolerance << ", max error " << error << std::endl;
        check(error < 500 * tolerance + 2e-8, "numeric trajectory agrees with analytic reduction");
        check(error < previous * .1 + 2e-9,
              "tightening integration must not worsen the comparison");
        previous = error;
    }
}
void periodIdentity() {
    const double amplitude = .6, e = std::cos(amplitude), delta = std::sqrt(2.);
    const AnalyticState x{amplitude, .5, 0, .5};
    const double h = invariants(config(), x).first;
    const double time =
        pi * std::sqrt((delta * delta - e * e) / (h * (1 + e))) *
        special::lauricellaFD(.5, {-.5, .5, -.5}, 1,
                              {(1 - e) / (delta - e), (e - 1) / (1 + e), (e - 1) / (delta + e)});
    const auto orbit = createSystem(config())->analyticSolutions()[0].instantiate(x, 0);
    const auto half = orbit.evaluate(time), whole = orbit.evaluate(2 * time);
    near(half[0] - half[2], -amplitude, 2e-10, "Lauricella relative half-period");
    near(whole[0] - whole[2], amplitude, 2e-10, "Lauricella relative full-period");
    for (double gravity : {0., 9.8}) {
        auto c = config();
        c.parameters["gravity_m_s2"] = gravity;
        const auto system = createSystem(c);
        const auto knowledge = system->dynamicalKnowledge();
        check(knowledge.invariantGuides.size() == (gravity == 0 ? 2 : 1),
              "Invariant applicability follows gravity");
        const AnalyticState y{.7, .2, -.3, .6};
        AnalyticState field(4);
        system->rhs(0, y.data(), field.data());
        for (const auto &guide : knowledge.invariantGuides) {
            double derivative = 0;
            for (std::size_t i = 0; i < 4; ++i) {
                auto plus = y, minus = y;
                plus[i] += 1e-6;
                minus[i] -= 1e-6;
                derivative += (guide.value(plus) - guide.value(minus)) / 2e-6 * field[i];
            }
            near(derivative, 0, 2e-8, "System invariant gradient annihilates native field");
        }
    }
}

void functions() {
    near(special::lauricellaFD(1, {1}, 2, {.7}), -std::log(.3) / .7, 2e-12,
         "elementary FD reduction");
    for (double m : {-5., -.2, 0., .7, .99})
        near(special::lauricellaFD(.5, {.5}, 1, {m}),
             2 *
                 (m < 0 ? special::ellipticK(m / (m - 1)) / std::sqrt(1 - m)
                        : special::ellipticK(m)) /
                 pi,
             3e-12, "elliptic FD reduction");
    near(special::lauricellaFD(.5, {.2, .3}, 2, {1, 1}), 4 / pi, 3e-12,
         "combined endpoint beta integral");
    near(special::lauricellaFD(.7, {.2, -.3}, 1.8, {-.4, -.4}),
         special::lauricellaFD(.7, {-.1}, 1.8, {-.4}), 3e-12, "coincident Lauricella arguments");
    near(special::lauricellaFD(1, {}, 2, {}), 1, 0, "empty Lauricella product");
    bool rejected = false;
    try {
        special::lauricellaFD(.5, {.5}, 1, {1});
    } catch (const std::exception &) {
        rejected = true;
    }
    check(rejected, "divergent endpoint rejected");
    rejected = false;
    try {
        special::lauricellaFD(1, {1}, 2, {1.01});
    } catch (const std::exception &) {
        rejected = true;
    }
    check(rejected, "complex branch rejected explicitly");
}
} // namespace
int main() {
    try {
        functions();
        periodIdentity();
        for (const auto &initial : std::vector<AnalyticState>{{.5, 0, .7, 0},
                                                              {0, 1, 0, 1},
                                                              {pi, -.4, 0, -.4},
                                                              {.6, .5, 0, .5},
                                                              {.2, 1, 0, -1},
                                                              {.2, -1, 0, 1},
                                                              {2 * pi + .6, .5, 0, .5},
                                                              {0, 0, 0, 1},
                                                              {0, 0, 0, -1}})
            compare(config(), initial);
        compare(config(), {0, 1, 0, -1}, 8);
        compare(config(), {0, -1, 0, 1}, 8);
        compare(config(), {0, 1, 0, -1.000001}, 12);
        compare(config(), {0, 1, 0, -.999999}, 12);
        // Primary paper examples use relative second angle and velocity; convert to native absolute
        // angles.
        for (double psi : {70 * pi / 180, 88 * pi / 180, 160 * pi / 180})
            compare(config(.5, 1, .5, 1), {pi / 8, 1, pi / 8 + psi, 3}, 10);
        auto c = config();
        c.parameters["gravity_m_s2"] = 1;
        check(!createSystem(c)->analyticSolutions()[0].instantiate({.2, 1, 0, 0}, 0),
              "gravity excluded");
        c = config();
        c.parameters["damping2_n_m_s"] = .01;
        check(!createSystem(c)->analyticSolutions()[0].instantiate({.2, 1, 0, 0}, 0),
              "damping excluded");
        c = config();
        c.driver.amplitude = .1;
        check(!createSystem(c)->analyticSolutions()[0].instantiate({.2, 1, 0, 0}, 0),
              "forcing excluded");
        std::cout << "Gravity-free double pendulum, conserved quantities, native residuals and "
                     "Lauricella identities passed\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
