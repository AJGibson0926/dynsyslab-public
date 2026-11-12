#pragma once

#include <dynsyslab/core/Integrator.hpp>
#include <dynsyslab/core/System.hpp>
#include <dynsyslab/workflow/Json.hpp>

#include <string>
#include <vector>

namespace dynsyslab::workflow::detail {

struct PrescribedControl {
    dynsyslab::core::ScalarControlLaw law;
    std::vector<double> activeWindow;
    std::string settingsCoordinates{"model"};
};

[[nodiscard]] PrescribedControl buildPrescribedControl(
    const std::string& kind, const json::Value& settings,
    const dynsyslab::core::ODESystem& system);

}  // namespace dynsyslab::workflow::detail
