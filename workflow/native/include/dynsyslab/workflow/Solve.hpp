#pragma once

#include <dynsyslab/workflow/Json.hpp>
#include <dynsyslab/core/Execution.hpp>

#include <filesystem>
#include <cstddef>
#include <string>
#include <vector>

namespace dynsyslab::workflow {

struct SolveResult {
    std::string label;
    std::string system;
    std::vector<std::filesystem::path> archives;
    std::filesystem::path hierarchyManifest;
    std::size_t trajectories{};
    std::size_t samplesPerTrajectory{};
    std::size_t acceptedSteps{};
    std::size_t rejectedSteps{};

    [[nodiscard]] std::wstring summary() const;
};

struct OdeInitialConditionSet {
    std::size_t trajectoryCount{};
    std::size_t stateDimension{};
    std::vector<std::string> stateLabels;
    std::vector<double> values;
};

// Stable capability list used by the CLI and desktop clients. Sampling
// algorithms remain in the workflow library rather than presentation code.
[[nodiscard]] std::vector<std::string> supportedOdeInitialConditionModes();

// Resolves exactly the deterministic ODE initial conditions that executeSolve
// will use, without integrating or writing archives. UI previews and scientific
// execution therefore share one generator and seed policy.
[[nodiscard]] OdeInitialConditionSet resolveOdeInitialConditions(
    const json::Value& job);

// Machine-readable, archive-free response used by `dynsyslab-cli generate-ics`.
[[nodiscard]] json::Value odeInitialConditionSetJson(
    const OdeInitialConditionSet& set);

[[nodiscard]] SolveResult executeSolve(
    const json::Value& job, const std::filesystem::path& archiveRoot,
    const std::filesystem::path& progressPath = {},
    const dynsyslab::core::ExecutionContext* execution = nullptr);

}  // namespace dynsyslab::workflow
