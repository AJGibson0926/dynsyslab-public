#include <algorithm>
#include <cmath>
#include <dynsyslab/core/PhaseSpace.hpp>
#include <iostream>
#include <stdexcept>
using namespace dynsyslab::core;
namespace {
constexpr double tau = 6.283185307179586476925286766559;
void check(bool ok, const char *text) {
    if (!ok)
        throw std::runtime_error(text);
}
void near(double a, double b, double tol, const char *text) {
    check(std::isfinite(a) && std::abs(a - b) < tol, text);
}
// An exact periodic translation of the curved-saddle oracle, with a time-dependent
// stable rate. This tests the forced variational equation and curved section geometry.
class ForcedSaddle final : public ODESystem {
    SystemDescriptor d_;
    SystemConfiguration c_;

  public:
    ForcedSaddle() {
        d_.identifier = "forced-saddle-test";
        d_.stateDimension = 2;
        c_.identifier = d_.identifier;
        c_.driver = {1., 1., .37};
    }
    const SystemDescriptor &descriptor() const override { return d_; }
    const SystemConfiguration &configuration() const override { return c_; }
    std::vector<double> initialState() const override { return {.12, -.08}; }
    void rhs(double t, const double *x, double *f) const override {
        const double angle = tau * t + .37, a = .7 + .2 * std::cos(angle), b = .4;
        const double u = x[0] - .12 * std::cos(angle), v = x[1] - .08 * std::sin(angle);
        f[0] = -.12 * tau * std::sin(angle) - a * u + .5 * (a + b) * v + (a + 2 * b) * v * v;
        f[1] = .08 * tau * std::cos(angle) + b * v;
    }
    void rhsWithInput(double, const double *x, double, double *f) const override {
        f[0] = -x[0];
        f[1] = x[1];
    }
    std::vector<double> jacobian(double, const double *) const override { return {-1, 0, 0, 1}; }
    std::vector<double> inputMatrix(const double *) const override { return {1, 1}; }
    bool invalid(double, const double *x) const override {
        return !std::isfinite(x[0]) || !std::isfinite(x[1]) ||
               std::abs(x[0]) + std::abs(x[1]) > 100;
    }
};
PhaseSpaceOptions options() {
    PhaseSpaceOptions o;
    o.lower = {-1., -1.};
    o.upper = {1., 1.};
    o.dynamics = "periodic-forcing";
    o.forcingSectionTime = .137;
    o.periodicSeeds = {{{.1, .1}, 0.}};
    o.outputStep = .025;
    o.maximumStep = .025;
    o.horizon = 4;
    o.transientTime = 0;
    o.vectorGridPerAxis = 3;
    o.manifoldIterations = 22;
    o.seedOffset = 1e-4;
    return o;
}
void exactSaddle() {
    ForcedSaddle system;
    auto o = options();
    auto result = analyzePhaseSpace(system, o);
    check(result.fixedPoints.empty(), "Forcing-section fixed point was mislabeled as equilibrium");
    check(result.periodicOrbits.size() == 1, "Forced periodic response missing");
    const auto &orbit = result.periodicOrbits[0];
    check(orbit.stability == "saddle_periodic_response" && orbit.stableDimension == 1 &&
              orbit.unstableDimension == 1,
          "Forced Floquet classification incorrect");
    near(result.options.forcingPeriod, 1., 1e-14, "Native forcing period");
    near(orbit.period, 1., 1e-14, "Response period");
    double lo = 10, hi = 0;
    for (auto m : orbit.floquetMultipliers) {
        lo = std::min(lo, m.real());
        hi = std::max(hi, m.real());
    }
    near(lo, std::exp(-.7), 2e-7, "Stable forced multiplier");
    near(hi, std::exp(.4), 2e-7, "Unstable forced multiplier");
    for (std::size_t i = 0; i < orbit.curve.time.size(); ++i) {
        const double phase = tau * orbit.curve.time[i] + .37;
        near(orbit.curve.values[2 * i], .12 * std::cos(phase), 2e-7,
             "Periodic orbit x matches absolute forcing time");
        near(orbit.curve.values[2 * i + 1], .08 * std::sin(phase), 2e-7,
             "Periodic orbit y matches absolute forcing time");
    }
    near(orbit.curve.time.front(), o.forcingSectionTime, 1e-14, "Section origin was reset to zero");
    near(orbit.curve.time.back(), o.forcingSectionTime + 1., 1e-13, "Exact return schedule");
    const double x0 = .12 * std::cos(tau * o.forcingSectionTime + .37),
                 y0 = .08 * std::sin(tau * o.forcingSectionTime + .37);
    check(result.manifolds.size() == 4, "Four forced saddle section branches missing");
    for (const auto &curve : result.manifolds) {
        if (curve.status == "unresolved")
            std::cerr << curve.reason << " error=" << curve.refinementError << '\n';
        check(curve.status != "unresolved" && curve.values.size() > 8,
              "Section curve failed independent refinement");
        check(curve.time.empty() && curve.baseKind == "periodic_response",
              "Section curve given a misleading time parameter");
        for (std::size_t i = 0; i < curve.values.size(); i += 2) {
            const double u = curve.values[i] - x0, v = curve.values[i + 1] - y0;
            if (curve.kind == "stable_section_manifold")
                near(v, 0., 2e-5, "Stable forcing-section manifold oracle");
            else
                near(u, .5 * v + v * v, 2e-3, "Nonlinear unstable forcing-section manifold oracle");
        }
    }
    for (std::size_t i = 0; i < result.vectorField.points.size(); i += 2) {
        double f[2];
        system.rhs(o.forcingSectionTime, &result.vectorField.points[i], f);
        near(result.vectorField.vectors[i], f[0], 1e-12, "Vector field uses forcing phase");
    }
    o.manifolds = false;
    o.forcingSectionTime += 2.;
    o.forcingPeriods = 2;
    auto shifted = analyzePhaseSpace(system, o);
    near(shifted.periodicOrbits.at(0).period, 1., 1e-13,
         "Subharmonic request retained a false doubled period");
    near(shifted.periodicOrbits[0].sectionState[0], orbit.sectionState[0], 2e-7,
         "Equivalent section phase changes the return map");
    o.forcingPeriod = .9;
    bool rejected = false;
    try {
        validatePhaseSpaceOptions(system, o);
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    check(rejected, "Inconsistent forcing period was accepted");
}
void sharedClocks() {
    SystemConfiguration c;
    c.identifier = "reddy_szeri";
    c.parameters = {
        {"kappa", 1.}, {"liquid_sound_speed_m_s", 1480.}, {"acoustic_pressure_amplitude", 10.}};
    auto system = createSystem(c);
    const double period = phaseSpaceForcingPeriod(*system);
    const double frequency = system->configuration().parameters.at("drive_frequency_hz");
    near(period * system->coordinates().secondsPerModelTime, 1 / frequency, 1e-15,
         "Internal acoustic clock units");
    c.driver = {1., 2 * frequency, 0.};
    near(phaseSpaceForcingPeriod(*createSystem(c)), period, 1e-12,
         "Commensurate internal/external clocks");
    c.driver.frequencyHz = std::sqrt(2.) * frequency;
    bool rejected = false;
    try {
        (void)phaseSpaceForcingPeriod(*createSystem(c));
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    check(rejected, "Incommensurate clocks were assigned a false common forcing period");
}
void rpeSaddle() {
    SystemConfiguration c;
    c.identifier = "rpe";
    c.parameters = {{"vapor_pressure_pa", 105000.}};
    c.driver = {20., 100000., .31};
    auto system = createSystem(c);
    auto o = options();
    const double T = phaseSpaceForcingPeriod(*system);
    o.lower = {.8, -.5};
    o.upper = {4., .5};
    o.forcingSectionTime = 0;
    o.periodicSeeds = {{{2.8, 0.}, T}};
    o.horizon = 2 * T;
    o.outputStep = T / 128;
    o.maximumStep = T / 64;
    o.manifoldIterations = 12;
    const auto result = analyzePhaseSpace(*system, o);
    check(result.periodicOrbits.size() == 1, "Tensile RPE saddle response missing");
    const auto &orbit = result.periodicOrbits[0];
    check(orbit.stability == "saddle_periodic_response" && result.manifolds.size() == 4,
          "Tensile RPE saddle branches missing");
    double integral = 0;
    for (std::size_t i = 1; i < orbit.curve.time.size(); ++i) {
        const double a = orbit.curve.values[2 * (i - 1)], b = orbit.curve.values[2 * i];
        integral +=
            (orbit.curve.time[i] - orbit.curve.time[i - 1]) * .5 * (1 / (a * a) + 1 / (b * b));
    }
    // Liouville: trace(Df)=-viscosity/r^2-3*v/r. The latter integral is an endpoint log.
    const double initial = orbit.curve.values[0],
                 finalRadius = orbit.curve.values[orbit.curve.values.size() - 2];
    const double determinant =
        std::exp(-system->configuration().dimensionlessGroups.at("viscosity") * integral -
                 3 * std::log(finalRadius / initial));
    near((orbit.floquetMultipliers[0] * orbit.floquetMultipliers[1]).real(), determinant, 2e-6,
         "RPE Floquet determinant violates the independent Liouville identity");
    for (const auto &branch : result.manifolds)
        check(branch.status != "unresolved" && branch.values.size() > 20 &&
                  branch.refinementError < 1e-4,
              "RPE section manifold failed refinement");
}
void rpeClockAndResponse() {
    SystemConfiguration c;
    c.identifier = "rpe";
    c.driver = {20., 100000., .31};
    auto system = createSystem(c);
    auto o = options();
    o.lower = {.8, -.2};
    o.upper = {1.2, .2};
    const double T = phaseSpaceForcingPeriod(*system);
    near(T * system->coordinates().secondsPerModelTime, 1. / 100000., 1e-15,
         "RPE physical/model forcing clock");
    o.forcingPeriod = T;
    o.forcingSectionTime = .17 * T;
    o.periodicSeeds = {{{1., 0.}, T}};
    o.horizon = 2 * T;
    o.maximumStep = T / 64;
    o.outputStep = T / 128;
    o.manifolds = false;
    auto result = analyzePhaseSpace(*system, o);
    if (result.periodicOrbits.empty())
        for (auto &s : result.diagnostics)
            std::cerr << s << '\n';
    check(result.periodicOrbits.size() == 1,
          "Forced RPE small-amplitude periodic response missing");
    check(result.periodicOrbits[0].closureResidual < 2e-7, "Forced RPE return closure");
    check(result.periodicOrbits[0].stability == "attracting_periodic_response",
          "Viscous forced RPE stability");
    double motion = 0;
    for (std::size_t i = 0; i < result.periodicOrbits[0].curve.values.size(); i += 2)
        motion = std::max(motion, std::abs(result.periodicOrbits[0].curve.values[i] - 1.));
    check(motion > 1e-7, "RPE forcing was silently ignored");
    ExecutionContext execution;
    execution.cancellation = []() { return true; };
    bool cancelled = false;
    try {
        (void)analyzePhaseSpace(*system, o, &execution);
    } catch (const ExecutionCancelled &) {
        cancelled = true;
    }
    check(cancelled, "Forced sidecar ignored cancellation");
}
} // namespace
int main() {
    try {
        exactSaddle();
        sharedClocks();
        rpeClockAndResponse();
        rpeSaddle();
        std::cout << "Periodic sidecar tests: PASS\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
