#pragma once

#include "dynsyslab_viewer/Theme.hpp"

#include <windows.h>

#include <optional>
#include <string>

namespace dsl::compute {

// Edits the high-use, controller-specific settings while preserving uncommon
// advanced keys. Q/R, KROM, representation, clock, and MPC policy remain on
// the primary Solver panel.
std::optional<std::string> showControllerSettingsDialog(
    HWND owner, const std::string& controllerKind,
    const std::string& currentSettings, const dsl::theme::Fonts& fonts);

}  // namespace dsl::compute
