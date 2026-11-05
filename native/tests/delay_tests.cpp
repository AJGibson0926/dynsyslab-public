#include <dynsyslab/core/Delay.hpp>

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
}

int main() {
    try {
        dynsyslab::core::DelayIntegrationOptions options;
        options.finalTime = 1.0;
        options.step = 0.05;
        const auto configuration = dynsyslab::core::resolveMackeyGlass();
        const auto result = dynsyslab::core::integrateMackeyGlass(
            configuration, {1.2, 1.0}, options);
        require(result.trajectoryCount == 2 && result.sampleCount == 21,
                "delay dimensions");
        require(result.time.front() == 0.0 && std::abs(result.time.back() - 1.0) < 1e-14,
                "delay time grid");
        const double expected = 1.2 + 0.05
            * (0.2 * 1.2 / (1.0 + std::pow(1.2, 10.0)) - 0.1 * 1.2);
        require(std::abs(result.values[1] - expected) < 1e-14,
                "Mackey-Glass Euler parity");
        require(dynsyslab::core::mackeyGlassDescriptor().family == "dde",
                "delay descriptor family");
        const auto catalog = dynsyslab::core::availableSystems();
        bool found = false;
        for (const auto& item : catalog) found = found || item.identifier == "mackey_glass";
        require(found, "delay catalog registration");
        dynsyslab::core::ExecutionContext limited;
        limited.limits.maximumMemoryBytes = 32;
        bool rejected = false;
        try {
            static_cast<void>(dynsyslab::core::integrateMackeyGlass(
                configuration, {1.2}, options, &limited));
        } catch (const dynsyslab::core::ExecutionLimitExceeded&) {
            rejected = true;
        }
        require(rejected, "delay memory limit");
        std::cout << "dynsyslab_delay_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_delay_tests: FAIL: " << error.what() << '\n';
        return 1;
    }
}
