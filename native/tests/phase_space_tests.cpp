#include <cmath>
#include <dynsyslab/core/PhaseSpace.hpp>
#include <iostream>
#include <stdexcept>
using namespace dynsyslab::core;
namespace {
void check(bool value, const char *message) {
    if (!value)
        throw std::runtime_error(message);
}
void near(double a, double b, double tolerance, const char *message) {
    check(std::isfinite(a) && std::abs(a - b) < tolerance, message);
}
// Grothendieck Labs CODE.dynamical-systems.dynamical-systems-model-library:
// Teschl's curved saddle and the class's exact radial-cycle oracle.
class Fixture final : public ODESystem {
    SystemDescriptor d_;
    SystemConfiguration c_;
    bool cycle_, curved_;

  public:
    Fixture(bool cycle, bool curved = false) : cycle_(cycle), curved_(curved) {
        d_.identifier = "phase-test";
        d_.stateDimension = 2;
        d_.preferredCoordinates = "dimensionless";
        c_.identifier = d_.identifier;
    }
    const SystemDescriptor &descriptor() const override { return d_; }
    const SystemConfiguration &configuration() const override { return c_; }
    std::vector<double> initialState() const override { return {.7, .1}; }
    DynamicalKnowledge dynamicalKnowledge() const override {
        DynamicalKnowledge knowledge;
        if (!cycle_ && !curved_) {
            knowledge.equilibriumGuides = [](const auto &, const auto &) {
                return std::vector<EquilibriumGuide>{{"deliberately-wrong-equilibrium", {.5, .5}}};
            };
            knowledge.invariantGuides.push_back({"deliberately-wrong-invariant",
                                                 "x is not conserved",
                                                 [](const auto &x) { return x[0]; }});
        }
        return knowledge;
    }
    void rhs(double t, const double *x, double *f) const override { rhsWithInput(t, x, 0, f); }
    void rhsWithInput(double, const double *x, double, double *f) const override {
        if (cycle_) {
            double radial = 1 - x[0] * x[0] - x[1] * x[1];
            f[0] = radial * x[0] - x[1];
            f[1] = x[0] + radial * x[1];
        } else if (curved_) {
            f[0] = -x[0] + x[1] + 3 * x[1] * x[1];
            f[1] = x[1];
        } else {
            f[0] = x[0] - x[0] * x[0] * x[0];
            f[1] = -x[1];
        }
    }
    std::vector<double> jacobian(double, const double *x) const override {
        if (cycle_)
            return {1 - 3 * x[0] * x[0] - x[1] * x[1], -1 - 2 * x[0] * x[1], 1 - 2 * x[0] * x[1],
                    1 - x[0] * x[0] - 3 * x[1] * x[1]};
        if (curved_)
            return {-1, 1 + 6 * x[1], 0, 1};
        return {1 - 3 * x[0] * x[0], 0, 0, -1};
    }
    std::vector<double> inputMatrix(const double *) const override { return {0, 0}; }
    bool invalid(double, const double *x) const override {
        return !std::isfinite(x[0]) || !std::isfinite(x[1]) || std::abs(x[0]) > 1e6 ||
               std::abs(x[1]) > 1e6;
    }
};
PhaseSpaceOptions options() {
    PhaseSpaceOptions o;
    o.lower = {-2, -2};
    o.upper = {2, 2};
    o.horizon = 12;
    o.transientTime = 10;
    o.periodicOrbits = false;
    return o;
}
void fixedAndManifolds() {
    Fixture system(false);
    auto o = options();
    auto result = analyzePhaseSpace(system, o);
    check(result.fixedPoints.size() == 3, "three known equilibria");
    check(result.manifolds.size() == 4, "four saddle branches");
    for (const auto &p : result.fixedPoints) {
        near(p.state[1], 0, 1e-8, "equilibria y=0");
        near(p.residual, 0, 1e-8, "equilibrium residual");
    }
    for (const auto &branch : result.manifolds) {
        check(branch.time.size() > 10, "nonempty manifold");
        check(branch.refinementError < 1e-3, "manifold seed-halving check");
        for (std::size_t i = 0; i < branch.time.size(); ++i)
            near(branch.values[i * 2 + (branch.kind == "stable_manifold" ? 0 : 1)], 0, 1e-9,
                 "exact separatrix axes");
    }
    check(result.vectorField.points.size() == 242, "sparse 11 x 11 vector grid");
    bool seedOnly = false, discrepancy = false;
    for (const auto &validation : result.knowledgeChecks) {
        if (validation.factId == "deliberately-wrong-equilibrium")
            seedOnly = validation.status == "seed-only";
        if (validation.factId == "deliberately-wrong-invariant" &&
            validation.status == "discrepancy")
            discrepancy = true;
    }
    check(seedOnly, "Incorrect researched equilibrium was taken for granted");
    check(discrepancy, "Incorrect researched invariant was taken for granted");
    Fixture curved(false, true);
    result = analyzePhaseSpace(curved, o);
    check(result.fixedPoints.size() == 1 && result.manifolds.size() == 4, "curved saddle branches");
    for (const auto &branch : result.manifolds)
        if (branch.kind == "unstable_manifold")
            for (std::size_t i = 0; i < branch.time.size(); ++i) {
                double x = branch.values[2 * i], y = branch.values[2 * i + 1];
                near(x, y / 2 + y * y, 2e-7, "curved manifold invariance oracle");
            }
}
void cycles() {
    Fixture system(true);
    auto o = options();
    o.manifolds = false;
    o.periodicOrbits = true;
    o.periodicSeeds = {{{0, -1.05}, 6.1}};
    o.sectionSeeds = {{0, -1}};
    auto result = analyzePhaseSpace(system, o);
    if (result.periodicOrbits.empty())
        for (const auto &d : result.diagnostics)
            std::cerr << d << '\n';
    check(result.periodicOrbits.size() == 1, "exact radial limit cycle found");
    auto orbit = result.periodicOrbits.front();
    near(orbit.period, 2 * 3.14159265358979323846, 1e-6, "exact cycle period");
    check(orbit.isolated && orbit.stability == "attracting_limit_cycle", "Floquet stability");
    check(orbit.sectionMultipliers.size() == 1,
          "Autonomous cycle did not expose its reduced section spectrum");
    near(orbit.sectionMultipliers.front().real(), std::exp(-4 * 3.14159265358979323846),
         1e-8, "Projected Poincare multiplier");
    check(orbit.transversalityAngle > .99 && orbit.tangentIdentityError < 1e-8,
          "Poincare transversality or neutral tangent identity failed");
    check(result.poincareMaps.size() == 1 && result.poincareMaps[0].status == "complete",
          "Explicit first-return map missing");
    near(result.poincareMaps[0].returnTime, 2 * 3.14159265358979323846, 2e-6,
         "First positive return time");
    near(result.poincareMaps[0].multipliers.front().real(),
         std::exp(-4 * 3.14159265358979323846), 2e-8,
         "First-return derivative projection");
    auto affine = o;
    affine.periodicOrbits = false;
    affine.fixedPoints = false;
    affine.vectorField = false;
    affine.sectionNormal = {1, 1};
    affine.sectionReference = {0, 0};
    const double inverseRootTwo = 1 / std::sqrt(2.0);
    affine.sectionSeeds = {{inverseRootTwo, -inverseRootTwo}};
    const auto affineResult = analyzePhaseSpace(system, affine);
    check(affineResult.poincareMaps.size() == 1 &&
              affineResult.poincareMaps[0].status == "complete",
          "Affine first-return map missing");
    near(affineResult.poincareMaps[0].returnTime, 2 * 3.14159265358979323846, 2e-6,
         "Affine first positive return time");
    near(affineResult.poincareMaps[0].multipliers.front().real(),
         std::exp(-4 * 3.14159265358979323846), 2e-8,
         "Affine first-return derivative projection");
    double transverse = 1;
    for (auto m : orbit.floquetMultipliers)
        if (std::abs(m - 1.) > .1)
            transverse = m.real();
    near(transverse, std::exp(-4 * 3.14159265358979323846), 1e-8, "exact transverse multiplier");
    for (std::size_t i = 0; i < orbit.curve.time.size(); ++i)
        near(std::hypot(orbit.curve.values[2 * i], orbit.curve.values[2 * i + 1]), 1, 2e-7,
             "exact invariant circle");
    o.periodicSeeds.clear();
    o.sectionSeeds.clear();
    o.gridPerAxis = 3;
    result = analyzePhaseSpace(system, o);
    check(!result.periodicOrbits.empty(), "oriented return search discovers cycle");
    SystemConfiguration c;
    c.identifier = "harmonic_oscillator";
    auto harmonic = createSystem(c);
    o.periodicSeeds = {{{0, 1}, 2 * 3.14159265358979323846}};
    result = analyzePhaseSpace(*harmonic, o);
    check(result.periodicOrbits.size() == 1 && !result.periodicOrbits[0].isolated,
          "center orbit is not an isolated limit cycle");
}
void shellHistory() {
    SystemConfiguration c;
    c.identifier = "marmottant";
    auto system = createSystem(c);
    auto o = options();
    o.lower = {.8, -.1};
    o.upper = {1.5, .1};
    bool rejected = false;
    try {
        validatePhaseSpaceOptions(*system, o);
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    check(rejected, "Shell history was projected into a two-state vector field");
    o.lower = {.8, -.1};
    o.upper = {.95, .1};
    validatePhaseSpaceOptions(*system, o);
}
void bubbleBranches() {
    SystemConfiguration c;c.identifier="marmottant";
    c.parameters={{"critical_surface_tension_n_m",.01}};
    auto system=createSystem(c);auto o=options();o.lower={1.01,-.1};o.upper={1.05,.1};
    bool rejected=false;
    try{validatePhaseSpaceOptions(*system,o);}catch(const std::invalid_argument&){rejected=true;}
    check(rejected,"Ruptured-shell guide crossed the later elastic boundary");
    c={};c.identifier="reddy_szeri";system=createSystem(c);
    o.lower=PhaseState(12,-.1);o.upper=PhaseState(12,.1);o.lower[0]=.8;o.upper[0]=.99;o.vectorField=false;
    validatePhaseSpaceOptions(*system,o);
    o.upper[0]=1.2;rejected=false;
    try{validatePhaseSpaceOptions(*system,o);}catch(const std::invalid_argument&){rejected=true;}
    check(rejected,"Gas-exponent switching was treated as a globally smooth field");
    c.parameters={{"kappa",1},{"liquid_sound_speed_m_s",1480},{"acoustic_pressure_amplitude",10}};
    system=createSystem(c);rejected=false;
    try{validatePhaseSpaceOptions(*system,o);}catch(const std::invalid_argument&){rejected=true;}
    check(rejected,"Internal acoustic translation forcing was treated as autonomous");
    c.parameters["translation_enabled"]=0;validatePhaseSpaceOptions(*createSystem(c),o);
    c.parameters["translation_enabled"]=1;c.parameters["drive_frequency_hz"]=0;
    validatePhaseSpaceOptions(*createSystem(c),o);
}
void cancellation() {
    Fixture system(false);
    auto o = options();
    ExecutionContext execution;
    execution.cancellation = []() { return true; };
    bool cancelled = false;
    try {
        (void)analyzePhaseSpace(system, o, &execution);
    } catch (const ExecutionCancelled &) {
        cancelled = true;
    }
    check(cancelled, "phase analysis cancellation");
}
} // namespace
int main() {
    try {
        fixedAndManifolds();
        cycles();
        cancellation();
        shellHistory();
        bubbleBranches();
        std::cout << "phase-space tests: PASS\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
