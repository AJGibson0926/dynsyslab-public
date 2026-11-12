#pragma once

#include <dynsyslab/core/Execution.hpp>
#include <dynsyslab/core/System.hpp>

#include <cstddef>
#include <vector>

namespace dynsyslab::core {

struct DelayIntegrationOptions {
    double initialTime{};
    double finalTime{100.0};
    double step{0.05};
};

struct DelaySolution {
    std::vector<double> time;
    std::vector<double> values;
    std::size_t trajectoryCount{};
    std::size_t sampleCount{};
    std::size_t acceptedSteps{};
};

// x' = -gamma*x + beta*x(t-delay), with constant prehistory; finite delayed-exponential sum.
[[nodiscard]] double linearDelayConstantHistory(double beta, double gamma, double delay,
                                                double history, double elapsed);
[[nodiscard]] std::vector<AnalyticSolutionFamily> mackeyGlassAnalyticSolutions(
    const SystemConfiguration& configuration);
[[nodiscard]] DynamicalKnowledge mackeyGlassKnowledge();
[[nodiscard]] const SystemDescriptor& mackeyGlassDescriptor();
[[nodiscard]] SystemConfiguration resolveMackeyGlass(
    const SystemConfiguration& requested = {});
[[nodiscard]] DelaySolution integrateMackeyGlass(
    const SystemConfiguration& configuration,
    const std::vector<double>& historyValues,
    const DelayIntegrationOptions& options,
    const ExecutionContext* execution = nullptr);

}  // namespace dynsyslab::core
