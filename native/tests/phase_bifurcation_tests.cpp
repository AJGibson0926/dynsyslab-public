#include <cmath>
#include <dynsyslab/core/PhaseSpace.hpp>
#include <iostream>
#include <memory>
#include <stdexcept>

using namespace dynsyslab::core;
namespace {
enum class NormalForm { fold, transcritical, pitchfork };

class ScalarNormalForm final : public ODESystem {
  public:
    ScalarNormalForm(NormalForm form, double parameter) : form_(form) {
        descriptor_.identifier = "scalar-normal-form-test";
        descriptor_.stateDimension = 1;
        configuration_.identifier = descriptor_.identifier;
        configuration_.parameters["mu"] = parameter;
    }
    const SystemDescriptor &descriptor() const override { return descriptor_; }
    const SystemConfiguration &configuration() const override { return configuration_; }
    std::vector<double> initialState() const override { return {0}; }
    void rhs(double, const double *x, double *f) const override {
        const double mu = configuration_.parameters.at("mu");
        if (form_ == NormalForm::fold)
            f[0] = mu - x[0] * x[0];
        else if (form_ == NormalForm::transcritical)
            f[0] = mu * x[0] - x[0] * x[0];
        else
            f[0] = mu * x[0] - x[0] * x[0] * x[0];
    }
    void rhsWithInput(double t, const double *x, double, double *f) const override {
        rhs(t, x, f);
    }
    std::vector<double> jacobian(double, const double *x) const override {
        const double mu = configuration_.parameters.at("mu");
        if (form_ == NormalForm::fold)
            return {-2 * x[0]};
        if (form_ == NormalForm::transcritical)
            return {mu - 2 * x[0]};
        return {mu - 3 * x[0] * x[0]};
    }
    std::vector<double> inputMatrix(const double *) const override { return {0}; }
    bool invalid(double, const double *x) const override { return !std::isfinite(x[0]); }

  private:
    NormalForm form_;
    SystemDescriptor descriptor_;
    SystemConfiguration configuration_;
};

class ForcedResponse final : public ODESystem {
  public:
    explicit ForcedResponse(double amplitude) : amplitude_(amplitude) {
        descriptor_.identifier = "forced-response-test";
        descriptor_.stateDimension = 1;
        configuration_.identifier = descriptor_.identifier;
        configuration_.driver = {1, 1, 0};
    }
    const SystemDescriptor &descriptor() const override { return descriptor_; }
    const SystemConfiguration &configuration() const override { return configuration_; }
    std::vector<double> initialState() const override { return {amplitude_}; }
    void rhs(double time, const double *x, double *f) const override {
        constexpr double tau = 6.283185307179586476925286766559;
        f[0] = -(x[0] - amplitude_ * std::cos(tau * time)) -
               amplitude_ * tau * std::sin(tau * time);
    }
    void rhsWithInput(double, const double *x, double, double *f) const override {
        f[0] = -x[0];
    }
    std::vector<double> jacobian(double, const double *) const override { return {-1}; }
    std::vector<double> inputMatrix(const double *) const override { return {0}; }
    bool invalid(double, const double *x) const override { return !std::isfinite(x[0]); }

  private:
    double amplitude_;
    SystemDescriptor descriptor_;
    SystemConfiguration configuration_;
};

void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}

PhaseBifurcationBranch run(NormalForm form, double initial, double step,
                           double seed, std::size_t maximumPoints = 40) {
    PhaseSpaceOptions options;
    options.lower = {-2};
    options.upper = {2};
    options.vectorField = false;
    options.fixedPoints = false;
    options.manifolds = false;
    options.periodicOrbits = false;
    options.residualTolerance = 1e-10;
    PhaseBifurcationOptions bifurcation;
    bifurcation.enabled = true;
    bifurcation.branch = "equilibrium";
    bifurcation.parameter = "mu";
    bifurcation.minimum = -.6;
    bifurcation.maximum = .8;
    bifurcation.initial = initial;
    bifurcation.step = step;
    bifurcation.minimumStep = 1e-4;
    bifurcation.maximumStep = .08;
    bifurcation.maximumPoints = maximumPoints;
    bifurcation.seedState = {seed};
    PhaseSystemFactory factory = [form](double mu) {
        return std::make_unique<ScalarNormalForm>(form, mu);
    };
    return continuePhaseBifurcation(factory, options, bifurcation);
}

bool has(const PhaseBifurcationBranch &branch, const std::string &kind) {
    for (const auto &event : branch.events)
        if (event.kind == kind)
            return true;
    return false;
}

void exactNormalForms() {
    const auto transcritical = run(NormalForm::transcritical, -.45, .05, 0, 28);
    require(transcritical.points.size() > 10 && has(transcritical, "transcritical"),
            "Transcritical structural classifier failed");
    const auto pitchfork = run(NormalForm::pitchfork, -.45, .05, 0, 28);
    require(pitchfork.points.size() > 10 && has(pitchfork, "pitchfork"),
            "Reflection-symmetric pitchfork classifier failed");
    const auto fold = run(NormalForm::fold, .5, -.04, std::sqrt(.5), 55);
    require(fold.points.size() > 15 && has(fold, "fold"),
            "Pseudo-arclength continuation did not cross and classify the fold");
    bool turned = false;
    for (std::size_t j = 1; j < fold.points.size(); ++j)
        if (fold.points[j].parameter > fold.points[j - 1].parameter)
            turned = true;
    require(turned, "Fold branch did not turn in the physical parameter");
}

void periodicBranches() {
    constexpr double tau = 6.283185307179586476925286766559;
    PhaseSpaceOptions autonomous;
    autonomous.lower = {-2, -2};
    autonomous.upper = {2, 2};
    autonomous.vectorField = false;
    autonomous.manifolds = false;
    autonomous.outputStep = .02;
    autonomous.maximumStep = .02;
    autonomous.horizon = 8;
    autonomous.sectionAxis = 0;
    autonomous.sectionValue = 0;
    autonomous.sectionDirection = 1;
    PhaseBifurcationOptions orbit;
    orbit.enabled = true;
    orbit.branch = "periodic-orbit";
    orbit.parameter = "mu";
    orbit.minimum = .6;
    orbit.maximum = 1.4;
    orbit.initial = 1;
    orbit.step = .04;
    orbit.minimumStep = 1e-4;
    orbit.maximumStep = .08;
    orbit.maximumPoints = 8;
    orbit.seedState = {0, -1};
    orbit.seedPeriod = tau;
    PhaseSystemFactory radial = [](double mu) {
        SystemConfiguration configuration;
        configuration.identifier = "radial_limit_cycle";
        configuration.parameters = {{"mu", mu}, {"alpha", 1}, {"omega", 1}};
        return createSystem(configuration);
    };
    const auto cycles = continuePhaseBifurcation(radial, autonomous, orbit);
    require(cycles.points.size() == orbit.maximumPoints,
            "Autonomous periodic-orbit continuation stopped early");
    for (const auto &point : cycles.points) {
        require(std::abs(std::hypot(point.state[0], point.state[1]) -
                         std::sqrt(point.parameter)) < 2e-6,
                "Radial periodic branch missed the exact cycle family");
        require(point.spectrum.size() == 1 &&
                    std::abs(point.spectrum[0].real() - std::exp(-2 * point.parameter * tau)) < 2e-7,
                "Autonomous branch used the neutral multiplier instead of the section spectrum");
    }

    PhaseSpaceOptions forced;
    forced.lower = {-2};
    forced.upper = {2};
    forced.dynamics = "periodic-forcing";
    forced.vectorField = false;
    forced.manifolds = false;
    forced.horizon = 4;
    forced.outputStep = 1. / 128;
    forced.maximumStep = 1. / 64;
    orbit.minimum = .2;
    orbit.maximum = 1.2;
    orbit.initial = .5;
    orbit.step = .05;
    orbit.seedState = {.5};
    orbit.seedPeriod = 0;
    orbit.maximumPoints = 6;
    PhaseSystemFactory response = [](double amplitude) {
        return std::make_unique<ForcedResponse>(amplitude);
    };
    const auto responses = continuePhaseBifurcation(response, forced, orbit);
    require(responses.points.size() == orbit.maximumPoints,
            "Forced periodic-response continuation stopped early");
    for (const auto &point : responses.points) {
        require(std::abs(point.state[0] - point.parameter) < 2e-7 &&
                    point.spectrum.size() == 1 &&
                    std::abs(point.spectrum[0].real() - std::exp(-1.)) < 2e-7,
                "Fixed-forcing-phase branch lost its exact response or multiplier");
    }
}
} // namespace

int main() {
    try {
        exactNormalForms();
        periodicBranches();
        std::cout << "phase bifurcation tests: PASS\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
