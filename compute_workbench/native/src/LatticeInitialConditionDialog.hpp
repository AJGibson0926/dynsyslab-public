#pragma once

#include "dynsyslab_viewer/Theme.hpp"

#include <windows.h>

#include <optional>
#include <string>

namespace dsl::compute {

// Generates or validates explicit Ising/Potts site arrays without expanding
// tens of thousands of sites into individual desktop controls.
std::optional<std::string> showLatticeInitialConditionDialog(
    HWND owner, const std::string& system, std::size_t size,
    std::size_t states, std::size_t trajectoryCount,
    const std::string& current, const dsl::theme::Fonts& fonts);

}  // namespace dsl::compute
