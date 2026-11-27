#include <dynsyslab/core/Integrator.hpp>

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

dynsyslab::core::NativeSolution solve(dynsyslab::core::IntegrationMethod method) {
    dynsyslab::core::SystemConfiguration configuration;
    configuration.identifier = "duffing";
    configuration.preset = "double_well_reference";
    auto system = dynsyslab::core::createSystem(configuration);
    dynsyslab::core::IntegrationOptions options;
    options.finalTime = 0.1;
    options.outputStep = 0.001;
    options.method = method;
    return dynsyslab::core::integrate(*system, {0.1, 0.0}, 1, options);
}

void testFixedAndAdaptive() {
    const auto fixed = solve(dynsyslab::core::IntegrationMethod::rk4);
    const auto adaptive = solve(dynsyslab::core::IntegrationMethod::dormandPrince45);
    require(fixed.sampleCount == 101 && adaptive.sampleCount == 101, "integration sample count");
    require(std::abs(fixed.state(0, 100)[0] - adaptive.state(0, 100)[0]) < 1.0e-10,
            "fixed/adaptive Duffing parity");
    require(adaptive.statistics.rightHandSideCalls > 0, "adaptive statistics missing");
    const auto series = fixed.trajectorySeries(0);
    require(series.valuesData() == fixed.values.data() && series.sampleCount() == fixed.sampleCount,
            "integrator trajectory TimeSeries must be a zero-copy view");
    require(!series.fourier().cached(), "integrator ran Fourier analysis eagerly");
}

void testStiffBDF() {
#ifdef DYNSYSLAB_SUNDIALS
    const auto reference = solve(dynsyslab::core::IntegrationMethod::dormandPrince45);
    const auto stiff = solve(dynsyslab::core::IntegrationMethod::bdf);
    require(stiff.sampleCount == 101, "BDF sample count");
    require(std::abs(reference.state(0, 100)[0] - stiff.state(0, 100)[0]) < 1.0e-7,
            "BDF/RK45 Duffing parity");
    require(stiff.statistics.acceptedSteps > 0 && stiff.statistics.rightHandSideCalls > 0,
            "BDF statistics missing");
#endif
}

void testBatchControlAndProgress() {
    dynsyslab::core::SystemConfiguration configuration;
    configuration.identifier = "van_der_pol";
    auto system = dynsyslab::core::createSystem(configuration);
    dynsyslab::core::IntegrationOptions options;
    options.finalTime = 0.02;
    options.outputStep = 0.001;
    options.method = dynsyslab::core::IntegrationMethod::rk4;
    double finalProgress = 0.0;
    auto solution = dynsyslab::core::integrate(*system, {0.1, 0.0, 0.2, 0.0}, 2, options,
        [](std::size_t trajectory, double, const double*) { return trajectory ? -0.1 : 0.1; },
        [&](const auto& update) { finalProgress = update.fraction; });
    require(solution.trajectoryCount == 2 && solution.controlValues.size() == 42,
            "controlled batch shape");
    require(solution.control(0, 0) == 0.1 && solution.control(1, 0) == -0.1,
            "control history values");
    require(finalProgress == 1.0, "progress did not complete");
}

void testDimensionlessRPE() {
    dynsyslab::core::SystemConfiguration configuration;
    configuration.identifier = "rpe";
    configuration.preset = "air_water_20c";
    auto system = dynsyslab::core::createSystem(configuration);
    dynsyslab::core::IntegrationOptions options;
    options.finalTime = 0.01;
    options.outputStep = 0.001;
    auto solution = dynsyslab::core::integrate(*system, {1.0, 0.0}, 1, options);
    require(solution.coordinateSystem == "dimensionless", "RPE coordinate metadata");
    require(std::abs(solution.state(0, 10)[0] - 1.0) < 1.0e-12, "RPE equilibrium drift");
}

void testSampledControlClock(dynsyslab::core::IntegrationMethod method) {
    dynsyslab::core::SystemConfiguration configuration;
    configuration.identifier = "van_der_pol";
    auto system = dynsyslab::core::createSystem(configuration);
    dynsyslab::core::IntegrationOptions options;
    options.finalTime = 0.3;
    options.outputStep = 0.05;
    options.initialStep = 0.01;
    options.method = method;
    std::size_t calls = 0;
    dynsyslab::core::ControlPolicy control;
    control.evaluation = dynsyslab::core::ControlEvaluationMode::sampled;
    control.sampleTime = 0.1;
    control.value = [&calls](std::size_t, double, const double*) {
        return static_cast<double>(++calls);
    };
    const auto solution = dynsyslab::core::integrateControlled(
        *system, {0.1, 0.0}, 1, options, std::move(control));
    require(calls == 4 && solution.statistics.controllerCalls == 4,
            "sampled controller runs once per clock tick");
    const std::vector<double> expected{1, 1, 2, 2, 3, 3, 4};
    require(solution.controlValues.size() == expected.size(),
            "sampled controller output shape");
    for (std::size_t index = 0; index < expected.size(); ++index)
        require(solution.controlValues[index] == expected[index],
                "sampled controller is right-continuous and held");
}

void testTruncateOnFailure(dynsyslab::core::IntegrationMethod method) {
    dynsyslab::core::SystemConfiguration configuration;
    configuration.identifier = "van_der_pol";
    auto system = dynsyslab::core::createSystem(configuration);
    dynsyslab::core::IntegrationOptions options;
    options.finalTime = 0.5; options.outputStep = 0.05; options.initialStep = 0.01;
    options.method = method;
    const auto policy = [] {
        dynsyslab::core::ControlPolicy control;
        control.evaluation = dynsyslab::core::ControlEvaluationMode::sampled;
        control.sampleTime = 0.05;
        // Trajectory 0 fails at t = 0.2 (a controller blowing up); trajectory 1 runs through.
        control.value = [](std::size_t trajectory, double time, const double*) {
            if (trajectory == 0 && time > 0.19) throw std::runtime_error("controller produced a nonfinite input");
            return 0.0;
        };
        return control;
    };
    bool aborted = false;
    try { (void)dynsyslab::core::integrateControlled(*system, {0.1, 0.0, 0.2, 0.0}, 2, options, policy()); }
    catch (const std::runtime_error&) { aborted = true; }
    require(aborted, "by default a failure aborts the solve");
    options.truncateOnFailure = true;
    const auto solution = dynsyslab::core::integrateControlled(*system, {0.1, 0.0, 0.2, 0.0}, 2, options, policy());
    require(solution.completedSamples.size() == 2 && solution.completedSamples[1] == solution.sampleCount,
            "the second trajectory completes");
    const auto kept = solution.completedSamples[0];
    require(kept >= 3 && kept < solution.sampleCount, "the failed trajectory keeps its computed samples");
    require(std::isfinite(solution.state(0, kept - 1)[0]) && std::isnan(solution.state(0, kept)[0]),
            "samples after the failure are NaN");
    require(solution.terminationReasons[0].find("nonfinite") != std::string::npos
            && solution.terminationReasons[1].empty(), "the failure reason is recorded");
}

void testInvalidSampledControlClock() {
    dynsyslab::core::SystemConfiguration configuration;
    configuration.identifier = "van_der_pol";
    auto system = dynsyslab::core::createSystem(configuration);
    dynsyslab::core::IntegrationOptions options;
    dynsyslab::core::ControlPolicy control;
    control.value = [](std::size_t, double, const double*) { return 0.0; };
    control.evaluation = dynsyslab::core::ControlEvaluationMode::sampled;
    bool rejected = false;
    try {
        (void)dynsyslab::core::integrateControlled(
            *system, {0.1, 0.0}, 1, options, std::move(control));
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    require(rejected, "sampled controller requires a positive sample time");
}

void testCancellation() {
    dynsyslab::core::SystemConfiguration configuration;
    configuration.identifier = "duffing";
    configuration.preset = "double_well_reference";
    auto system = dynsyslab::core::createSystem(configuration);
    dynsyslab::core::IntegrationOptions options;
    options.finalTime = 10.0;
    options.outputStep = 0.001;
    dynsyslab::core::ExecutionContext execution;
    std::size_t probes = 0;
    execution.cancellation = [&] { return ++probes > 10; };
    bool cancelled = false;
    try {
        static_cast<void>(dynsyslab::core::integrate(
            *system, system->initialState(), 1, options, {}, {}, &execution));
    } catch (const dynsyslab::core::ExecutionCancelled&) {
        cancelled = true;
    }
    require(cancelled, "integration cancellation checkpoint");
}

void testMemoryLimit() {
    dynsyslab::core::SystemConfiguration configuration;
    configuration.identifier = "duffing";
    auto system = dynsyslab::core::createSystem(configuration);
    dynsyslab::core::IntegrationOptions options;
    options.finalTime = 1.0;
    options.outputStep = 0.001;
    dynsyslab::core::ExecutionContext execution;
    execution.limits.maximumMemoryBytes = 1024;
    bool limited = false;
    try {
        static_cast<void>(dynsyslab::core::integrate(
            *system, system->initialState(), 1, options, {}, {}, &execution));
    } catch (const dynsyslab::core::ExecutionLimitExceeded&) {
        limited = true;
    }
    require(limited, "integration memory limit");
}

}  // namespace

int main() {
    try {
        testFixedAndAdaptive();
        testStiffBDF();
        testBatchControlAndProgress();
        testDimensionlessRPE();
        testSampledControlClock(dynsyslab::core::IntegrationMethod::rk4);
#ifdef DYNSYSLAB_SUNDIALS
        testSampledControlClock(dynsyslab::core::IntegrationMethod::bdf);
#endif
        testTruncateOnFailure(dynsyslab::core::IntegrationMethod::rk4);
#ifdef DYNSYSLAB_SUNDIALS
        testTruncateOnFailure(dynsyslab::core::IntegrationMethod::bdf);
#endif
        testInvalidSampledControlClock();
        testCancellation();
        testMemoryLimit();
        std::cout << "dynsyslab_integrator_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_integrator_tests: FAIL: " << error.what() << '\n';
        return 1;
    }
}
