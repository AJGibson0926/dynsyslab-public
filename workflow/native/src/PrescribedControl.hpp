#pragma once

#include <dynsyslab/core/Integrator.hpp>
#include <dynsyslab/core/System.hpp>
#include <dynsyslab/workflow/Json.hpp>

#include <functional>
#include <string>
#include <vector>

namespace dynsyslab::workflow::detail {

struct PrescribedControl {
    dynsyslab::core::ScalarControlLaw law;
    std::vector<double> activeWindow;
    std::string settingsCoordinates{"model"};
    // Optional per-trajectory design record (tape digest, pair_id/polarity, ...).
    std::function<json::Value(std::size_t)> trajectoryMetadata;
};

// waveform, random_pulse, tabulated, command_tape, telegraph_multisine, symmetric_pulse, pulse_branches.
[[nodiscard]] bool isPrescribedControl(const std::string& kind);

[[nodiscard]] PrescribedControl buildPrescribedControl(
    const std::string& kind, const json::Value& settings,
    const dynsyslab::core::ODESystem& system);

// The designed-input kinds of PrescribedInputs.cpp (all but waveform and random_pulse).
[[nodiscard]] PrescribedControl buildDesignedInput(
    const std::string& kind, const json::Value& settings,
    const dynsyslab::core::CoordinateScaling& scaling);

}  // namespace dynsyslab::workflow::detail
