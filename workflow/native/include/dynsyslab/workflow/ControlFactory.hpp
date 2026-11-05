#pragma once

#include <dynsyslab/core/Integrator.hpp>
#include <dynsyslab/core/LQR.hpp>
#include <dynsyslab/core/System.hpp>
#include <dynsyslab/workflow/Json.hpp>

#include <cstddef>
#include <string>
#include <vector>

namespace dynsyslab::workflow {

// A resolved controller plus the target and provenance needed by the archive.
// State-space LQR/SDRE and Koopman KLQR remain separate controller families.
struct ControllerSetup {
    dynsyslab::core::ScalarControlPolicy policy;
    dynsyslab::core::StateTargetFunction target; // optional state reference for plotting/provenance
    dynsyslab::core::StateTargetFunction observableTarget;
    std::vector<std::string> targetObservableNames;
    std::string targetSpace{"state"};
    std::string kind{"none"};
    std::string family{"none"};
    std::string evaluation{"none"};
    std::string settingsCoordinates{"model"};
    double sampleTime{};
    json::Value modelBindingEvidence;

    explicit operator bool() const noexcept {
        return static_cast<bool>(policy);
    }
};

[[nodiscard]] ControllerSetup buildController(
    const json::Value& controllerSpec, const json::Value& solverSpec,
    const dynsyslab::core::ODESystem& system, std::size_t trajectoryCount,
    double defaultSampleTime, const std::vector<double>& linearizationState);

}  // namespace dynsyslab::workflow
