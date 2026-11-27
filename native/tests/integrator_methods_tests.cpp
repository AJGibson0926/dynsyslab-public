// Native ODE methods and solve features of the phase-B port: RK23, DOP853, Radau IIA, CVODE
// Adams, LSODA switching, arbitrary output times, guards, step budget, checkpoints and pause,
// additive actuation, dimensional integration, PCHIP/nearest resampling, delay rollout and the
// adaptive method of lines.

#include <dynsyslab/core/Delay.hpp>
#include <dynsyslab/core/DelayRollout.hpp>
#include <dynsyslab/core/Lattice.hpp>
#include <dynsyslab/core/DimensionalSystem.hpp>
#include <dynsyslab/core/Integrator.hpp>
#include <dynsyslab/core/PDE.hpp>
#include <dynsyslab/core/TimeSeriesOperations.hpp>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace core = dynsyslab::core;

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

std::unique_ptr<core::ODESystem> system(const char* identifier, const char* preset = "") {
    core::SystemConfiguration configuration;
    configuration.identifier = identifier;
    configuration.preset = preset;
    return core::createSystem(configuration);
}

// Harmonic oscillator x'' = -x from (1, 0): x(t) = cos t.
core::NativeSolution harmonic(core::IntegrationMethod method, double rtol = 1.0e-10) {
    core::SystemConfiguration configuration;
    configuration.identifier = "harmonic_oscillator";
    configuration.parameters["damping_kg_s"] = 0.0;
    auto model = core::createSystem(configuration);
    core::IntegrationOptions options;
    options.finalTime = 10.0;
    options.outputStep = 0.5;
    options.method = method;
    options.relativeTolerance = rtol;
    options.absoluteTolerance = rtol * 1.0e-2;
    const auto initial = model->initialState();
    std::vector<double> start(initial.size(), 0.0);
    start[0] = 1.0;
    return core::integrate(*model, start, 1, options);
}

void testMethodsAgree() {
    const auto reference = harmonic(core::IntegrationMethod::dormandPrince853, 1.0e-12);
    const double exact = reference.state(0, reference.sampleCount - 1)[0];
    for (const auto method : {core::IntegrationMethod::bogackiShampine23, core::IntegrationMethod::dormandPrince45,
                              core::IntegrationMethod::radau5}) {
        const auto solution = harmonic(method, method == core::IntegrationMethod::bogackiShampine23 ? 1.0e-9 : 1.0e-10);
        const double value = solution.state(0, solution.sampleCount - 1)[0];
        require(std::abs(value - exact) < 2.0e-6, core::integrationMethodName(method) + " disagrees with DOP853");
    }
    // Undamped unit oscillator from (1, 0): x = cos t, and x^2 + v^2 stays 1.
    require(std::abs(exact - std::cos(10.0)) < 1.0e-9, "DOP853 accuracy");
    const auto* last = reference.state(0, reference.sampleCount - 1);
    require(std::abs(last[0] * last[0] + last[1] * last[1] - 1.0) < 1.0e-8, "DOP853 energy drift");
    // DOP853 needs far fewer steps than RK23 at a tight tolerance.
    const auto low = harmonic(core::IntegrationMethod::bogackiShampine23, 1.0e-9);
    require(reference.statistics.acceptedSteps < low.statistics.acceptedSteps, "DOP853 step economy");
    require(core::integrationMethodFromName("dop853") == core::IntegrationMethod::dormandPrince853
            && core::integrationMethodFromName("Radau") == core::IntegrationMethod::radau5
            && core::integrationMethodFromName("LSODA") == core::IntegrationMethod::lsoda,
            "method names");
}

// Van der Pol with mu = 1000 is stiff: Radau must finish within a small step budget.
void testStiffRadau() {
    core::SystemConfiguration configuration;
    configuration.identifier = "van_der_pol";
    configuration.parameters["mu"] = 1000.0;
    auto model = core::createSystem(configuration);
    core::IntegrationOptions options;
    options.finalTime = 1.0;
    options.outputStep = 0.25;
    options.method = core::IntegrationMethod::radau5;
    options.relativeTolerance = 1.0e-6;
    options.absoluteTolerance = 1.0e-8;
    options.maximumInternalSteps = 2000;
    const auto solution = core::integrate(*model, {2.0, 0.0}, 1, options);
    require(std::isfinite(solution.state(0, solution.sampleCount - 1)[0]), "stiff Radau finite");
    require(std::abs(solution.state(0, solution.sampleCount - 1)[0] - 2.0) < 0.05, "stiff Radau slow manifold");
#ifdef DYNSYSLAB_SUNDIALS
    options.method = core::IntegrationMethod::bdf;
    const auto bdf = core::integrate(*model, {2.0, 0.0}, 1, options);
    require(std::abs(bdf.state(0, bdf.sampleCount - 1)[0] - solution.state(0, solution.sampleCount - 1)[0]) < 1.0e-3,
            "Radau and BDF agree on the stiff problem");
    options.method = core::IntegrationMethod::lsoda;
    options.maximumInternalSteps = 100000;
    const auto switching = core::integrate(*model, {2.0, 0.0}, 1, options);
    require(std::abs(switching.state(0, switching.sampleCount - 1)[0] - solution.state(0, solution.sampleCount - 1)[0]) < 1.0e-3,
            "LSODA agrees on the stiff problem");
    const auto adams = harmonic(core::IntegrationMethod::adams, 1.0e-10);
    require(std::abs(adams.state(0, adams.sampleCount - 1)[0] - std::cos(10.0)) < 1.0e-6, "CVODE Adams accuracy");
    (void)harmonic(core::IntegrationMethod::bdf, 1.0e-8);
    const auto nonstiff = harmonic(core::IntegrationMethod::lsoda, 1.0e-10);
    require(nonstiff.statistics.methodSwitches == 0, "LSODA stays on Adams for a nonstiff problem");
#endif
}

void testOutputTimesGuardsBudget() {
    auto model = system("duffing", "double_well_reference");
    core::IntegrationOptions options;
    options.initialTime = 0.0; options.finalTime = 1.0;
    options.outputTimes = {0.0, 0.1, 0.35, 0.36, 1.0};
    options.method = core::IntegrationMethod::dormandPrince853;
    options.relativeTolerance = 1.0e-12; options.absoluteTolerance = 1.0e-14;
    const auto solution = core::integrate(*model, {0.1, 0.0}, 1, options);
    require(solution.sampleCount == 5 && solution.time[2] == 0.35, "arbitrary output times");
    core::IntegrationOptions uniform = options;
    uniform.outputTimes.clear(); uniform.outputStep = 0.05;
    const auto reference = core::integrate(*model, {0.1, 0.0}, 1, uniform);
    require(std::abs(reference.state(0, 20)[0] - solution.state(0, 4)[0]) < 1.0e-9, "output clock does not change the trajectory");

    auto guarded = uniform;
    guarded.method = core::IntegrationMethod::rk4;
    guarded.stateUpperBounds = {0.1000001, 1.0};
    guarded.stateLowerBounds = {-1.0, -1.0};
    guarded.truncateOnFailure = true;
    const auto clipped = core::integrate(*model, {0.1, 0.05}, 1, guarded);
    require(clipped.completedSamples[0] < clipped.sampleCount
            && clipped.terminationReasons[0].find("upper bound") != std::string::npos, "state bound guard");

    auto budget = uniform;
    budget.maximumInternalSteps = 3;
    budget.truncateOnFailure = true;
    const auto starved = core::integrate(*model, {0.1, 0.0, 0.1, 0.0}, 2, budget);
    require(starved.terminationReasons[0].find("step budget") != std::string::npos
            && starved.terminationReasons[1].find("step budget") != std::string::npos, "per-trajectory step budget");
}

void testCheckpointPause() {
    auto model = system("duffing", "double_well_reference");
    core::IntegrationOptions options;
    options.finalTime = 1.0; options.outputStep = 0.1;
    options.method = core::IntegrationMethod::rk4;
    options.checkpointInterval = 0.3;
    std::vector<double> seen;
    options.checkpoint = [&](const core::IntegrationCheckpoint& point) {
        seen.push_back(point.time);
        require(point.values[point.sample * point.stateDimension] == point.state[0], "checkpoint rows");
        return point.time < 0.55;   // pause at the second checkpoint
    };
    const auto solution = core::integrate(*model, {0.1, 0.0, 0.2, 0.0}, 2, options);
    require(solution.paused && seen.size() == 2, "pause at checkpoint");
    require(solution.completedSamples[0] == 7 && solution.completedSamples[1] == 0, "paused sample counts");
    require(solution.terminationReasons[1] == "paused", "unstarted trajectory marked paused");
}

void testAdditiveActuation() {
    auto model = system("duffing", "double_well_reference");
    core::IntegrationOptions options;
    options.finalTime = 0.5; options.outputStep = 0.1;
    options.method = core::IntegrationMethod::rk4;
    core::ControlPolicy additive;
    additive.vector = [](std::size_t, double, const double*, double* u) { u[0] = 0.3; u[1] = -0.2; };
    additive.inputDimension = 2;
    additive.additiveMatrix = {1.0, 0.0, 0.0, 1.0};
    const auto forced = core::integrateControlled(*model, {0.1, 0.0}, 1, options, additive);
    require(forced.inputDimension == 2 && forced.controlSeries(0).componentCount() == 2, "vector control history");
    const auto free = core::integrate(*model, {0.1, 0.0}, 1, options);
    require(std::abs(forced.state(0, 1)[0] - free.state(0, 1)[0]) > 1.0e-3, "additive input moves the state");
    core::ControlPolicy bad;
    bad.vector = additive.vector; bad.inputDimension = 2;
    bool rejected = false;
    try { (void)core::integrateControlled(*model, {0.1, 0.0}, 1, options, bad); }
    catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "vector input without actuation matrix is rejected");
}

void testDimensionalIntegration() {
    auto model = system("rpe");
    const auto scale = model->coordinates();
    require(core::supportsDimensionalIntegration(*model), "RPE declares physical scales");
    core::IntegrationOptions tau;
    tau.finalTime = 2.0; tau.outputStep = 0.5;
    tau.method = core::IntegrationMethod::dormandPrince853;
    tau.relativeTolerance = 1.0e-11; tau.absoluteTolerance = 1.0e-13;
    const std::vector<double> start{1.05, 0.0};
    const auto modelSolution = core::integrate(*model, start, 1, tau);
    auto physical = core::systemInCoordinates(system("rpe"), "dimensional");
    require(physical->descriptor().preferredCoordinates == "dimensional", "dimensional descriptor");
    const double T = scale.secondsPerModelTime, R0 = scale.physicalStatePerModelState[0];
    core::IntegrationOptions seconds = tau;
    seconds.finalTime = 2.0 * T; seconds.outputStep = 0.5 * T;
    seconds.absoluteTolerance = 1.0e-13 * R0;
    const auto dimensional = core::integrate(*physical, {1.05 * R0, 0.0}, 1, seconds);
    const double modelRadius = modelSolution.state(0, modelSolution.sampleCount - 1)[0];
    const double physicalRadius = dimensional.state(0, dimensional.sampleCount - 1)[0] / R0;
    require(std::abs(modelRadius - physicalRadius) < 1.0e-7, "dimensional and dimensionless RPE agree");
}

void testResampling() {
    const std::vector<double> time{0.0, 1.0, 2.0, 3.0, 4.0};
    const std::vector<double> values{0.0, 0.0, 1.0, 1.0, 1.0};  // monotone step
    const core::TimeSeries series(time, values, 1);
    std::vector<double> fine;
    for (int k = 0; k <= 40; ++k) fine.push_back(0.1 * k);
    const auto pchip = series.resample(fine, {core::Interpolation::pchip, core::Extrapolation::error});
    const auto cubic = series.resample(fine, {core::Interpolation::cubic, core::Extrapolation::error});
    bool cubicOvershoots = false;
    for (std::size_t k = 0; k < fine.size(); ++k) {
        require(pchip.value(k, 0) >= -1.0e-15 && pchip.value(k, 0) <= 1.0 + 1.0e-15, "PCHIP stays within the data");
        if (k && pchip.value(k, 0) + 1.0e-15 < pchip.value(k - 1, 0)) throw std::runtime_error("PCHIP is not monotone");
        cubicOvershoots = cubicOvershoots || cubic.value(k, 0) > 1.0 + 1.0e-9 || cubic.value(k, 0) < -1.0e-9;
    }
    require(cubicOvershoots, "Catmull-Rom cubic overshoots a step (PCHIP does not)");
    const auto nearest = series.resample({1.4, 1.5, 1.6}, {core::Interpolation::nearest, core::Extrapolation::error});
    require(nearest.value(0, 0) == 0.0 && nearest.value(1, 0) == 0.0 && nearest.value(2, 0) == 1.0, "nearest ties go down");
    require(core::interpolationFromName("pchip") == core::Interpolation::pchip, "interpolation names");
    // PCHIP reproduces data exactly at the knots and a straight line exactly.
    const core::TimeSeries line(time, {0.0, 2.0, 4.0, 6.0, 8.0}, 1);
    const auto onLine = line.resample({0.25, 3.75}, {core::Interpolation::pchip, core::Extrapolation::error});
    require(std::abs(onLine.value(0, 0) - 0.5) < 1.0e-14 && std::abs(onLine.value(1, 0) - 7.5) < 1.0e-14, "PCHIP is exact on lines");
}

void testDelayRollout() {
    // x'(t) = -x(t - 1), x = 1 on [-1, 0]: x(t) = 1 - t on [0, 1] (method of steps).
    core::DelayRolloutOptions options;
    options.step = 0.01;
    options.lags = {0.0, 1.0};
    options.outputSamples = 101;
    const std::vector<double> history(101, 1.0);
    const auto rollout = core::rolloutDelayModel(
        [](const double* window, double* derivative) { derivative[0] = -window[1]; }, 1, history, options);
    require(rollout.samples == 101 && rollout.status == core::DelayRolloutStatus::complete, "delay rollout length");
    require(std::abs(rollout.values[100] - 0.0) < 1.0e-12 && std::abs(rollout.values[50] - 0.5) < 1.0e-12,
            "first method-of-steps interval is exact");
    options.outputSamples = 201;
    const auto second = core::rolloutDelayModel(
        [](const double* window, double* derivative) { derivative[0] = -window[1]; }, 1, history, options);
    // On [1, 2]: x = 1 - t + (t - 1)^2 / 2.
    require(std::abs(second.values[200] - (1.0 - 2.0 + 0.5)) < 1.0e-6, "second method-of-steps interval");
    options.magnitudeLimits = {0.75};
    const auto stopped = core::rolloutDelayModel(
        [](const double* window, double* derivative) { derivative[0] = window[1]; }, 1, history, options);
    require(stopped.status == core::DelayRolloutStatus::magnitudeGuard && stopped.samples < 201, "delay magnitude guard");
}

void testAdaptivePDE() {
    core::StructuredMesh mesh;
    mesh.shape = {16}; mesh.spacing = {1.0 / 16.0}; mesh.origin = {0.0};
    mesh.boundary = {core::PDEBoundary::periodic};
    std::vector<double> field(16);
    for (std::size_t i = 0; i < 16; ++i) field[i] = std::sin(6.283185307179586 * i / 16.0);
    const auto rhs = [](double, const core::StructuredMesh&, const double* u, std::size_t n, double* out) {
        const double h = 1.0 / static_cast<double>(n);
        for (std::size_t i = 0; i < n; ++i) out[i] = 0.01 * (u[(i + n - 1) % n] - 2.0 * u[i] + u[(i + 1) % n]) / (h * h);
    };
    core::PDEIntegrationOptions options;
    options.finalTime = 1.0; options.outputStep = 0.25;
    options.method = core::PDEIntegrationMethod::adaptive;
    options.adaptiveMethod = core::IntegrationMethod::dormandPrince853;
    options.relativeTolerance = 1.0e-10; options.absoluteTolerance = 1.0e-12;
    const auto adaptive = core::integratePDE(mesh, {16}, field, options, rhs);
    // Semi-discrete decay rate of the first mode: lambda = 0.01 * (2 cos(2 pi h) - 2) / h^2.
    const double h = 1.0 / 16.0, lambda = 0.01 * (2.0 * std::cos(6.283185307179586 * h) - 2.0) / (h * h);
    require(adaptive.sampleCount == 5, "adaptive PDE samples");
    require(std::abs(adaptive.values[4 * 16 + 4] - std::exp(lambda) * field[4]) < 1.0e-8, "adaptive MOL accuracy");
    options.adaptiveMethod = core::IntegrationMethod::radau5;
    const auto stiff = core::integratePDE(mesh, {16}, field, options, rhs);
    require(std::abs(stiff.values[4 * 16 + 4] - std::exp(lambda) * field[4]) < 1.0e-6, "Radau MOL accuracy");
}

// Durable mid-run history for delay and lattice solves.
void testDelayAndLatticeCheckpoints() {
    core::DelayIntegrationOptions options;
    options.finalTime = 10.0; options.step = 0.1; options.checkpointSteps = 25;
    std::vector<std::size_t> rows;
    double lastValue = 0.0;
    options.checkpoint = [&](std::size_t, std::size_t samples, const double* times, const double* values) {
        rows.push_back(samples); lastValue = values[samples - 1];
        require(std::abs(times[samples - 1] - 0.1 * static_cast<double>(samples - 1)) < 1.0e-12, "delay checkpoint clock");
    };
    const auto delay = core::integrateMackeyGlass({}, {1.2}, options);
    require(rows == std::vector<std::size_t>({26, 51, 76}), "delay checkpoints every 25 steps before the end");
    require(lastValue == delay.values[75], "delay checkpoint holds the computed samples");

    core::IsingRequest ising;
    ising.size = 6; ising.spins.assign(36, 1); ising.temperature = 2.5; ising.seed = 3; ising.sweeps = 10;
    std::vector<std::size_t> sweeps;
    ising.checkpoint.interval = 4;
    ising.checkpoint.callback = [&](std::size_t done, const std::vector<double>& observations,
                                    const std::vector<std::int64_t>& spins) {
        sweeps.push_back(done);
        require(observations.size() == done * 2 && spins.size() == 36, "lattice checkpoint shapes");
    };
    const auto run = core::runIsing2D(ising);
    require(sweeps == std::vector<std::size_t>({4, 8}) && run.observations.size() == 20, "lattice checkpoints");
}

}  // namespace

int main() {
    try {
        testMethodsAgree();
        testStiffRadau();
        testOutputTimesGuardsBudget();
        testCheckpointPause();
        testAdditiveActuation();
        testDimensionalIntegration();
        testResampling();
        testDelayRollout();
        testAdaptivePDE();
        testDelayAndLatticeCheckpoints();
        std::cout << "dynsyslab_integrator_methods_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_integrator_methods_tests: FAIL: " << error.what() << '\n';
        return 1;
    }
}
