#include <dynsyslab/core/Lattice.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

void requireNear(double value, double expected, const char* message) {
    if (std::abs(value - expected) > 1.0e-12) throw std::runtime_error(message);
}

void testIsing2D() {
    dynsyslab::core::IsingRequest request;
    request.size = 2;
    request.spins = {1, 1, 1, 1};
    request.temperature = 1.0e-9;
    request.seed = 7;
    request.sweeps = 2;
    const auto result = dynsyslab::core::runIsing2D(request);
    require(result.observations.size() == 4 && result.spins.size() == 4,
            "unexpected 2-D Ising result shape");
    require(result.initialSpins == request.spins, "2-D Ising initial snapshot");
    requireNear(result.observations[0], 1.0, "2-D Ising magnetization");
    requireNear(result.observations[1], -2.0, "2-D Ising energy density");
}

void testIsing3D() {
    dynsyslab::core::IsingRequest request;
    request.size = 2;
    request.dimensions = 3;
    request.spins.assign(8, -1);
    request.temperature = 1.0e-9;
    request.sweeps = 1;
    const auto result = dynsyslab::core::runIsing3D(request);
    require(result.spins.size() == 8, "unexpected 3-D Ising result shape");
    requireNear(result.observations[0], -1.0, "3-D Ising magnetization");
    requireNear(result.observations[1], -3.0, "3-D Ising energy density");
}

void testPotts() {
    dynsyslab::core::PottsRequest request;
    request.size = 2;
    request.states = 3;
    request.spins.assign(4, 0);
    request.temperature = 1.0e-9;
    request.sweeps = 1;
    const auto result = dynsyslab::core::runPotts2D(request);
    require(result.spins.size() == 4, "unexpected Potts result shape");
    requireNear(result.observations[0], 1.0, "Potts order parameter");
    requireNear(result.observations[1], -2.0, "Potts energy density");
}

void testRandomizationAndLimits() {
    dynsyslab::core::IsingRequest request;
    request.size = 4;
    request.temperature = 2.269;
    request.seed = 42;
    request.sweeps = 3;
    const auto first = dynsyslab::core::runIsing2D(request);
    const auto second = dynsyslab::core::runIsing2D(request);
    require(first.initialSpins == second.initialSpins
            && first.observations == second.observations,
            "native lattice seed reproducibility");
    dynsyslab::core::ExecutionContext limited;
    limited.limits.maximumMemoryBytes = 16;
    bool rejected = false;
    try {
        static_cast<void>(dynsyslab::core::runIsing2D(request, &limited));
    } catch (const dynsyslab::core::ExecutionLimitExceeded&) {
        rejected = true;
    }
    require(rejected, "native lattice memory limit");
}

void testCatalog() {
    dynsyslab::core::SystemConfiguration requested;
    requested.identifier = "ising_3d";
    const auto resolved = dynsyslab::core::resolveLatticeSystem(requested);
    require(resolved.identifier == "ising_3d"
            && resolved.parameters.at("temperature") == 4.5115,
            "3-D Ising defaults");
    const auto systems = dynsyslab::core::availableSystems();
    const auto found = std::find_if(systems.begin(), systems.end(), [](const auto& item) {
        return item.identifier == "potts" && item.family == "lattice";
    });
    require(found != systems.end(), "lattice system catalog registration");
}

}  // namespace

int main() {
    try {
        testIsing2D();
        testIsing3D();
        testPotts();
        testRandomizationAndLimits();
        testCatalog();
        std::cout << "dynsyslab_lattice_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_lattice_tests: FAIL: " << error.what() << '\n';
        return 1;
    }
}
