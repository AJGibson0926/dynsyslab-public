#include <dynsyslab/core/KROM.hpp>
#include <dynsyslab/core/Hierarchy.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/ODE.hpp>

#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

void testODE() {
    dynsyslab::core::FixedIntegrationRequest request;
    request.model = "duffing";
    request.parameters = {0.2, 1.0, 1.0, 0.0, 0.0, 0.0};
    request.initialStates = {0.1, 0.0};
    request.finalTime = 0.01;
    request.stepSize = 0.001;
    const auto result = dynsyslab::core::integrateFixed(request);
    require(result.samples == 11 && result.stateDimension == 2, "unexpected ODE result shape");
    require(std::isfinite(result.state(0, 10)[0]), "non-finite ODE result");
}

void testKROM() {
    using Complex = std::complex<double>;
    dynsyslab::core::KROMRolloutRequest request;
    request.eigenvalues = {Complex{-0.1, 1.0}};
    request.modes = {Complex{1.0, 0.0}, Complex{0.0, 1.0}};
    request.initialCoordinates = {Complex{1.0, 0.0}};
    request.time = {0.0, 0.5, 1.0};
    request.stateDimension = 2;
    const auto result = dynsyslab::core::rolloutKROM(request);
    require(result.lifted.size() == 3 && result.reconstructed.size() == 6,
            "unexpected KROM result shape");
}

void testHierarchy() {
    dynsyslab::core::SystemConfiguration configuration;
    configuration.identifier = "duffing";
    configuration.preset = "double_well_reference";
    configuration.parameters = {{"damping", 0.2}};
    const auto key = dynsyslab::core::parameterSetKey(configuration);
    require(key.size() == 16, "parameter-set key size");
    dynsyslab::core::BatchRecord batch;
    batch.id = "batch"; batch.label = "test"; batch.createdUtc = "now";
    dynsyslab::core::SystemRunRecord system;
    system.id = "system"; system.systemIdentifier = "duffing";
    system.parameterSetKey = key;
    dynsyslab::core::ClusterRecord cluster; cluster.id = "cluster";
    dynsyslab::core::TrajectoryRecord trajectory;
    trajectory.id = "trajectory"; trajectory.initialState = {0.1, 0.0};
    trajectory.solutionArchive = "trajectory.npz";
    trajectory.timeSeries.push_back({"state", "time", "values"});
    cluster.trajectories.push_back(std::move(trajectory));
    system.clusters.push_back(std::move(cluster)); batch.systems.push_back(std::move(system));
    const auto encoded = dynsyslab::core::hierarchyJson(batch);
    require(encoded.find("dynsyslab.hierarchy/v2") != std::string::npos,
            "hierarchy schema serialization");
    require(encoded.find("trajectory.npz") != std::string::npos,
            "hierarchy time-series reference");
}

void testJsonPaths() {
    auto value = dynsyslab::core::json::parse(
        R"({"controller":{"settings":{"frequencies":[1,2]}}})");
    dynsyslab::core::json::setPath(
        value, "controller.settings.frequencies.1",
        dynsyslab::core::json::Value(3.0));
    require(value.at("controller").at("settings").at("frequencies")
            .array()[1].number() == 3.0,
        "indexed JSON campaign path");
    bool rejected = false;
    try {
        dynsyslab::core::json::setPath(
            value, "controller.settings.frequencies.2",
            dynsyslab::core::json::Value(4.0));
    } catch (const std::out_of_range&) {
        rejected = true;
    }
    require(rejected, "out-of-range JSON campaign path rejection");
}

}  // namespace

int main() {
    try {
        testODE();
        testKROM();
        testHierarchy();
        testJsonPaths();
        std::cout << "dynsyslab_core_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_core_tests: FAIL: " << error.what() << '\n';
        return 1;
    }
}
