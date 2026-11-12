#include <algorithm>
#include <cmath>
#include <dynsyslab/core/Integrator.hpp>
#include <dynsyslab/core/System.hpp>
#include <iostream>
#include <stdexcept>
using namespace dynsyslab::core;
namespace {
void check(bool ok, const std::string &text) {
    if (!ok)
        throw std::runtime_error(text);
}
void near(double a, double b, double tolerance, const std::string &label) {
    if (!std::isfinite(a) || !std::isfinite(b) || std::abs(a - b) > tolerance * (1 + std::abs(b)))
        throw std::runtime_error(label + ": " + std::to_string(a) + " vs " + std::to_string(b));
}
AnalyticTrajectory orbit(const ODESystem &system, const std::string &id, const AnalyticState &x,
                         double t0 = 0) {
    for (const auto &f : system.analyticSolutions())
        if (f.descriptor.identifier == id)
            return f.instantiate(x, t0);
    throw std::runtime_error("Missing system-owned analytic family " + id);
}
void compare(const SystemConfiguration &c, const std::string &id, const AnalyticState &initial,
             double duration, double step = 1e-5) {
    auto system = createSystem(c);
    const double t0 = .371;
    const auto exact = orbit(*system, id, initial, t0);
    check(bool(exact), id + ": " + exact.reason);
    const auto start = exact.evaluate(t0);
    for (std::size_t j = 0; j < initial.size(); ++j)
        near(start[j], initial[j], 2e-13, "Initial condition");
    for (double f : {-.01, .01, .1, .3, .7, 1.}) {
        const double t = t0 + f * duration;
        const auto y = exact.evaluate(t), a = exact.evaluate(t - step),
                   b = exact.evaluate(t + step);
        AnalyticState field(y.size());
        system->resetRuntimeState();
        system->rhs(t, y.data(), field.data());
        for (std::size_t j = 0; j < y.size(); ++j)
            near((b[j] - a[j]) / (2 * step), field[j], 3e-6, id + " native RHS");
    }
    double previous = 1e100;
    for (double tolerance : {1e-7, 1e-10}) {
        IntegrationOptions options;
        options.initialTime = t0;
        options.finalTime = t0 + duration;
        options.outputStep = duration / 70;
        options.maximumStep = duration / 40;
        options.relativeTolerance = tolerance;
        options.absoluteTolerance = tolerance * .01;
        system->resetRuntimeState();
        const auto result = integrate(*system, initial, 1, options);
        double error = 0;
        for (std::size_t i = 0; i < result.time.size(); ++i) {
            const auto y = exact.evaluate(result.time[i]);
            for (std::size_t j = 0; j < y.size(); ++j)
                error = std::max(error, std::abs(y[j] - result.values[i * y.size() + j]));
        }
        std::cout << id << " tolerance " << tolerance << " error " << error << std::endl;
        check(error < 1500 * tolerance + 1e-8, id + " numerical agreement");
        check(error < .2 * previous + 2e-8, id + " numerical refinement");
        previous = error;
    }
}
void verifyGuides(const SystemConfiguration &c, const AnalyticState &lo, const AnalyticState &hi) {
    auto system = createSystem(c);
    const auto k = system->dynamicalKnowledge();
    check(k.researchStatus.find("reviewed") != std::string::npos, "Reviewed system metadata");
    check(!k.papers.empty(), "System owns bibliography");
    check(bool(k.equilibriumGuides), "Executable candidate guides");
    auto guides = k.equilibriumGuides(lo, hi);
    check(!guides.empty(), "At least one candidate");
    for (const auto &guide : guides) {
        AnalyticState field(lo.size());
        system->resetRuntimeState();
        system->rhs(0, guide.state.data(), field.data());
        for (double v : field)
            near(v, 0, 2e-7, "Candidate native residual");
    }
    for (const auto &invariant : k.invariantGuides) {
        AnalyticState x(lo.size()), field(lo.size());
        for (std::size_t j = 0; j < x.size(); ++j)
            x[j] = .37 * lo[j] + .63 * hi[j];
        system->rhs(0, x.data(), field.data());
        double defect = 0;
        for (std::size_t j = 0; j < x.size(); ++j) {
            auto a = x, b = x;
            a[j] -= 1e-6;
            b[j] += 1e-6;
            defect += (invariant.value(b) - invariant.value(a)) / 2e-6 * field[j];
        }
        near(defect, 0, 2e-6, "Native invariant directional derivative");
    }
}
void neural() {
    SystemConfiguration c;
    c.identifier = "hodgkin_huxley";
    verifyGuides(c, {-100, 0, 0, 0}, {50, 1, 1, 1});
    check(!createSystem(c)->dynamicalKnowledge().facts[1].applicable,
          "Ordinary free-voltage HH is not a constant-voltage gating limit");
    check(!orbit(*createSystem(c), "hh-constant-voltage-gates", {-65, .2, .3, .4}),
          "No artificial clamp of ordinary HH");
    for (double v : {-65., -40., -55.}) {
        c.parameters = {{"g_na", 5}, {"e_na", v}, {"g_k", 0}, {"g_l", 0}, {"bias_current", 0}};
        compare(c, "hh-constant-voltage-gates", {v, .2, .3, .4}, .003, 1e-8);
    }
    c.parameters["bias_current"] = 1;
    check(!orbit(*createSystem(c), "hh-constant-voltage-gates", {-55, .2, .3, .4}),
          "Current mismatch excludes gating reference");
}
void pairs() {
    SystemConfiguration c;
    c.identifier = "kuramoto";
    c.parameters = {{"oscillators", 2}, {"frequency_spread", .5}};
    for (double coupling : {-3., -1., -.2, 0., .4, 1., 3.}) {
        c.parameters["coupling"] = coupling;
        compare(c, "kuramoto-two-oscillators", {-.3, .4}, 20);
    }
    for (double spread : {-2., .001, 2.}) {
        c.parameters["frequency_spread"] = spread;
        c.parameters["coupling"] = .7;
        compare(c, "kuramoto-two-oscillators", {30, -8}, 10);
    }
    c.parameters["frequency_spread"] = .5;
    for (double coupling : {1 - 1e-6, 1 + 1e-6}) {
        c.parameters["coupling"] = coupling;
        compare(c, "kuramoto-two-oscillators", {0, 3}, 30);
    }
    c.parameters["oscillators"] = 1;
    compare(c, "kuramoto-uncoupled", {.7}, 2);
    c.parameters["oscillators"] = 2;
    c.parameters["coupling"] = 2;
    verifyGuides(c, {-4, -4}, {4, 4});
    c.parameters["oscillators"] = 37;
    check(createSystem(c)->descriptor().stateLabels.back() == R"(\theta_{37})",
          "Configured array labels extend past catalogue default");
    c.parameters["oscillators"] = 7;
    c.parameters["frequency_spread"] = 0;
    verifyGuides(c, AnalyticState(7, -1), AnalyticState(7, 1));
    c.parameters["oscillators"] = 2;
    c.parameters["coupling"] = -3;
    const auto stationary = orbit(*createSystem(c), "kuramoto-two-oscillators", {0, 0});
    check(stationary.evaluate(1e4) == AnalyticState({0, 0}),
          "Pure decaying projective eigenvector is preserved");
    c.parameters["frequency_spread"] = 2;
    c.parameters["coupling"] = .7;
    const auto system = createSystem(c);
    const auto x = AnalyticState{.2, -.4};
    const auto q = orbit(*system, "kuramoto-two-oscillators", x);
    const double period = 2 * 3.14159265358979323846 / std::sqrt(16 - .49);
    for (double t : {-100 * period, -.7, 0., .2, 20 * period}) {
        const auto y = q.evaluate(t), z = q.evaluate(t + period);
        near(z[0] - y[0], -3.14159265358979323846, 1e-11, "Phase-slip lift first phase");
        near(z[1] - y[1], 3.14159265358979323846, 1e-11, "Phase-slip lift second phase");
        const auto reverse = orbit(*system, "kuramoto-two-oscillators", y, t).evaluate(0);
        near(reverse[0], x[0], 1e-10, "Flow inverse");
        near(reverse[1], x[1], 1e-10, "Flow inverse");
    }
}
void springs() {
    SystemConfiguration c;
    c.identifier = "spring_coupled_pendulums";
    c.parameters = {{"damping1_n_m_s", 0},
                    {"damping2_n_m_s", 0},
                    {"spring_constant_n_m", 40},
                    {"mass1_kg", .7},
                    {"mass2_kg", 2}};
    verifyGuides(c, {-4, -1, -4, -1}, {4, 1, 4, 1});
    compare(c, "spring-pendulum-synchronous", {.7, 0, .7, 0}, 3);
    compare(c, "spring-pendulum-synchronous", {.1, 7, .1, 7}, 3);
    c.parameters["spring_rest_length_m"] = .8;
    check(!orbit(*createSystem(c), "spring-pendulum-synchronous", {.7, 0, .7, 0}),
          "Prestress excludes synchronous reference");
    check(createSystem(c)
              ->dynamicalKnowledge()
              .equilibriumGuides({-.1, -1, -.1, -1}, {.1, 1, .1, 1})
              .empty(),
          "Prestretched vertical rods are not equilibrium seeds");
}
void henonHeiles() {
    SystemConfiguration c;c.identifier="henon_heiles";
    for(double lam:{-2.,.3,1.}){
        c.parameters={{"lam",lam}};
        for(double angle:{0.,2.09439510239319549,4.18879020478639098}){
            const double ex=std::sin(angle),ey=std::cos(angle);
            compare(c,"henon-heiles-straight-line",{ex*.3/lam,ey*.3/lam,ex*-.1/lam,ey*-.1/lam},25);
        }
    }
    c.parameters={{"lam",1.}};
    compare(c,"henon-heiles-straight-line",{0,-.5,0,0},8);
    compare(c,"henon-heiles-straight-line",{0,.999,0,0},8);
    compare(c,"henon-heiles-straight-line",{0,1e-5,0,-2e-5},12);
    auto system=createSystem(c);
    const auto knowledge=system->dynamicalKnowledge();
    check(knowledge.facts.size()>1 && knowledge.facts[0].id!=knowledge.facts[1].id,"System research facts need distinct identities");
    check(!orbit(*system,"henon-heiles-straight-line",{1e-5,.3,0,.1}),"Near-axis motion must not be guessed as exact");
    check(!orbit(*system,"henon-heiles-straight-line",{0,1.7,0,0}),"Escaping branch must not wrap through an elliptic pole");
    c.driver={1.,1.,0};check(!orbit(*createSystem(c),"henon-heiles-straight-line",{0,.3,0,0}),"Forced HH is not the autonomous reference");
}
void bubbles() {
    SystemConfiguration c;
    c.identifier = "marmottant";
    c.parameters = {{"buckling_radius_m", .975e-6 * .99}};
    verifyGuides(c, {.995, -.1}, {1.015, .1});
    compare(c, "marmottant-stationary", {1, 0}, 2);
    check(!orbit(*createSystem(c), "marmottant-stationary", {1, .01}),
          "Moving shell is not stationary");
    check(createSystem(c)->dynamicalKnowledge().equilibriumGuides({.9, -1}, {1.2, 1}).empty(),
          "No cross-history shell guide");
    c = {};
    c.identifier = "reddy_szeri";
    c.parameters = {{"vapor_pressure_pa", 0}};
    AnalyticState lo(12, -1), hi(12, 1);
    lo[0] = .8;
    hi[0] = 1.2;
    verifyGuides(c, lo, hi);
    AnalyticState initial{1, 0, .3, 0, .1, .02, .05, 0, -.03, .01, .02, -.01};
    compare(c, "reddy-szeri-stationary-radius-shapes", initial, 3);
    c.parameters["vapor_pressure_pa"] = 2330;
    const auto native = createSystem(c);
    const auto guides = native->dynamicalKnowledge().equilibriumGuides(lo, hi);
    check(!guides.empty(), "Non-unit equilibrium found");
    initial[0] = guides[0].state[0];
    check(std::abs(initial[0] - 1) > 1e-4, "Do not assume nominal radius is equilibrium");
    compare(c, "reddy-szeri-stationary-radius-shapes", initial, 2);
    initial[0] = 1;
    check(!orbit(*native, "reddy-szeri-stationary-radius-shapes", initial),
          "Reject wrong radial pressure balance");
    initial[0] = guides[0].state[0];
    c.parameters["liquid_sound_speed_m_s"] = 1480;
    c.parameters["acoustic_pressure_amplitude"] = 10;
    check(!orbit(*createSystem(c), "reddy-szeri-stationary-radius-shapes", initial),
          "Reject effective acoustic translation forcing");
}
} // namespace
int main() {
    try {
        neural();
        pairs();
        springs();
        bubbles();
        henonHeiles();
        std::cout << "Catalogue exact limits and guides passed\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
