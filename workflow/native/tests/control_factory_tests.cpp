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


void testTargetRampToggle() {
    // A ramp works on any target kind; null/false/{"enabled":false} switch it off.
    const auto plain = dynsyslab::workflow::resolveControlTarget(dynsyslab::core::json::parse(
        R"({"kind":"constant","state":[1.5,0]})"), 2);
    close(plain(0.0)[0], 1.5, 0.0, "unramped constant target");
    for (const char* off : {R"({"kind":"constant","state":[1.5,0],"ramp":null})",
                            R"({"kind":"constant","state":[1.5,0],"ramp":false})",
                            R"({"kind":"constant","state":[1.5,0],"ramp":{"enabled":false,"duration":2,"from":[1,0]}})"})
        close(dynsyslab::workflow::resolveControlTarget(dynsyslab::core::json::parse(off), 2)(0.0)[0],
              1.5, 0.0, "disabled ramp changed the target");
    const auto ramped = dynsyslab::workflow::resolveControlTarget(dynsyslab::core::json::parse(
        R"({"kind":"constant","state":[1.5,0],"ramp":{"duration":2,"from":[1,0]}})"), 2);
    close(ramped(0.0)[0], 1.0, 1.0e-14, "ramp start");
    close(ramped(1.0)[0], 1.25, 1.0e-14, "ramp midpoint");
    close(ramped(1.0)[1], 0.5 * 30.0 * 0.0625 / 2.0, 1.0e-14, "ramp velocity is the exact derivative");
    close(ramped(3.0)[0], 1.5, 1.0e-14, "ramp end");
    // Waveform ramp: the velocity stays the derivative of the ramped position.
    const auto wave = dynsyslab::workflow::resolveControlTarget(dynsyslab::core::json::parse(
        R"({"kind":"waveform","offset":[1.5,0],"amplitudes":[0.4],"angular_frequencies":[1.9],
            "phases_rad":[0],"ramp":{"duration":10,"from":[1,0]}})"), 2);
    const double t = 3.7, h = 1.0e-5;
    close(wave(t)[1], (wave(t + h)[0] - wave(t - h)[0]) / (2 * h), 1.0e-7, "ramped waveform velocity");
    bool rejected = false;
    try { (void)dynsyslab::workflow::resolveControlTarget(dynsyslab::core::json::parse(
        R"({"kind":"constant","state":[1,0],"ramp":{"duration":0,"from":[1,0]}})"), 2); }
    catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "zero-length ramp accepted");
}

void testControlCapToggle() {
    using dynsyslab::workflow::controlCap;
    const auto parse = [](const char* text) { return dynsyslab::core::json::parse(text); };
    require(!controlCap(parse("null")) && !controlCap(parse("false")), "null/false must leave the controller uncapped");
    const auto symmetric = controlCap(parse("0.5"));
    require(symmetric && symmetric->first == -0.5 && symmetric->second == 0.5, "scalar cap");
    const auto interval = controlCap(parse("[-1.2,0.3]"));
    require(interval && interval->first == -1.2 && interval->second == 0.3, "asymmetric cap");
    for (const char* bad : {"0", "-1", "[1,0]", "[1,2,3]", "true"}) {
        bool rejected = false;
        try { (void)controlCap(parse(bad)); } catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "invalid control_limit accepted");
    }
    // The cap applies to classical controllers too, and is absent unless asked for.
    const auto system = rpe();
    const double state[]{5.0, -5.0};
    const auto uncapped = controller(R"({
      "kind":"lqr","settings":{"A":[0,1,-1,-1],"B":[0,1],"Q":100,"R":0.01},
      "target":{"kind":"constant","state":[0,0]}})", *system);
    const auto capped = controller(R"({
      "kind":"lqr","settings":{"A":[0,1,-1,-1],"B":[0,1],"Q":100,"R":0.01,"control_limit":[-0.25,0.1]},
      "target":{"kind":"constant","state":[0,0]}})", *system);
    const double free = uncapped.policy.value(0, 0.0, state), clipped = capped.policy.value(0, 0.0, state);
    require(std::abs(free) > 0.25, "test state does not exercise the cap");
    require(clipped == (free > 0 ? 0.1 : -0.25), "classical LQR cap not applied");
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
        testTargetRampToggle();
        testControlCapToggle();
        std::cout << "dynsyslab_control_factory_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_control_factory_tests: FAIL: " << error.what() << '\n';
        return 1;
    }
}
