#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/System.hpp>
#include <dynsyslab/workflow/ControlFactory.hpp>
#include <dynsyslab/workflow/ControlTargets.hpp>

#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

void close(double actual, double expected, double tolerance, const char* message) {
    if (std::abs(actual - expected) > tolerance)
        throw std::runtime_error(std::string(message) + ": got "
            + std::to_string(actual) + ", expected " + std::to_string(expected));
}

std::unique_ptr<dynsyslab::core::ODESystem> rpe() {
    dynsyslab::core::SystemConfiguration configuration;
    configuration.identifier = "rpe";
    configuration.preset = "air_water_20c";
    return dynsyslab::core::createSystem(configuration);
}

dynsyslab::workflow::ControllerSetup controller(
    const std::string& specification, const dynsyslab::core::ODESystem& system) {
    const auto source = dynsyslab::core::json::parse(specification);
    const auto solver = dynsyslab::core::json::parse(
        R"({"control_evaluation":"automatic"})");
    return dynsyslab::workflow::buildController(source, solver, system, 1, 0.001,
                                                 system.controlledInitialState());
}

void testPhysicalWaveformScaling() {
    const auto system = rpe();
    const auto scaling = system->coordinates();
    require(scaling.physicalTimeAvailable, "RPE physical clock unavailable");
    require(scaling.physicalInputAvailable, "RPE physical input unavailable");
    close(scaling.physicalInputPerModelInput, 100000.0, 1.0e-9,
          "RPE pressure scale");
    const double frequency = 1.0 / scaling.secondsPerModelTime;
    const std::string source = std::string(R"({
      "kind":"waveform","settings":{
        "frequencies_hz":[)") + std::to_string(frequency) + R"(],
        "amplitudes":[100000],"phases_rad":[0],
        "active_time_window":[0,)" + std::to_string(0.5 * scaling.secondsPerModelTime)
        + R"(]},"target":{"kind":"constant","state":[1,0]}})";
    const auto result = controller(source, *system);
    const double state[]{1.0, 0.0};
    close(result.policy.value(0, 0.25, state), 1.0, 2.0e-6,
          "physical waveform conversion");
    close(result.policy.value(0, 0.75, state), 0.0, 1.0e-14,
          "physical active-window conversion");
    require(result.settingsCoordinates == "physical",
            "physical waveform provenance");
    require(!result.target, "open-loop waveform must not manufacture a target");
}

void testMultitoneAndModelCoordinates() {
    const auto system = rpe();
    const double state[]{1.0, 0.0};
    const auto result = controller(R"({
      "kind":"waveform","settings":{
        "frequencies":[1,2],"model_amplitudes":[2,0.5],
        "phases_rad":[1.5707963267948966,1.5707963267948966],"model_offset":0.25},
      "target":{"kind":"constant","state":[1,0]}})", *system);
    close(result.policy.value(0, 0.0, state), 2.75, 1.0e-14,
          "model-coordinate multitone");
    require(result.settingsCoordinates == "model", "model waveform provenance");
}

void testRandomPulseScaling() {
    const auto system = rpe();
    const auto scaling = system->coordinates();
    const std::string source = std::string(R"({
      "kind":"random_pulse","settings":{"amplitude":100000,"pulse_duration":)")
        + std::to_string(scaling.secondsPerModelTime)
        + R"(,"seed":17},"target":{"kind":"constant","state":[1,0]}})";
    const auto result = controller(source, *system);
    const double state[]{1.0, 0.0};
    const double first = result.policy.value(0, 0.1, state);
    close(result.policy.value(0, 0.9, state), first, 0.0,
          "physical pulse-duration conversion");
    require(std::abs(first) <= 1.0, "physical pulse-amplitude conversion");
    require(result.settingsCoordinates == "physical", "physical pulse provenance");
    require(!result.target, "open-loop random forcing must not manufacture a target");
}

void testWaveformTarget() {
    const auto specification = dynsyslab::core::json::parse(R"({
      "kind":"waveform","offset":[1,0],"frequencies":[0.5],
      "amplitudes":[0.2],"phases_rad":[0]})");
    const auto target = dynsyslab::workflow::resolveControlTarget(specification, 2);
    const auto state = target(0.0);
    close(state[0], 1.0, 1.0e-14, "waveform target position");
    close(state[1], 0.2 * 3.14159265358979323846, 1.0e-14,
          "waveform target kinematic derivative");
}

void testTrajectoryTarget() {
    const auto path = std::filesystem::temp_directory_path()
        / "dynsyslab-control-target-test.npz";
    dynsyslab::core::NpzRecord archive;
    archive.arrays["time"] = dynsyslab::core::NpyArray::fromDoubles(
        {3}, {0.0, 1.0, 2.0});
    archive.arrays["values"] = dynsyslab::core::NpyArray::fromDoubles(
        {3, 2}, {0.0, 0.0, 1.0, 2.0, 2.0, 4.0});
    dynsyslab::core::writeNpzAtomic(path, archive);
    const auto encoded = std::string("{\"kind\":\"trajectory\",\"source\":\"")
        + path.generic_u8string() + "\"}";
    const auto target = dynsyslab::workflow::resolveControlTarget(
        dynsyslab::core::json::parse(encoded), 2);
    const auto state = target(0.5);
    close(state[0], 0.5, 1.0e-14, "trajectory target interpolation");
    close(state[1], 1.0, 1.0e-14, "trajectory target second component");
    std::error_code ignored; std::filesystem::remove(path, ignored);
}

void testExplicitClassicalLqrMatrices() {
    const auto system = rpe();
    const auto result = controller(R"({
      "kind":"lqr","settings":{"A":[0,1,-1,-1],"B":[0,1],"Q":1,"R":1},
      "target":{"kind":"constant","state":[0,0]}})", *system);
    const double state[]{0.1, -0.2};
    require(std::isfinite(result.policy.value(0, 0.0, state)),
            "explicit classical LQR matrices were not accepted");
}

}  // namespace

int main() {
    try {
        testPhysicalWaveformScaling();
        testMultitoneAndModelCoordinates();
        testRandomPulseScaling();
        testWaveformTarget();
        testTrajectoryTarget();
        testExplicitClassicalLqrMatrices();
        std::cout << "dynsyslab_control_factory_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_control_factory_tests: FAIL: " << error.what() << '\n';
        return 1;
    }
}
