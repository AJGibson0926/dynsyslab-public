#pragma once

#include "dynsyslab_viewer/Theme.hpp"

#include <windows.h>

#include <optional>
#include <string>

namespace dsl::compute {

// Returns a canonical JSON object containing only explicit parameter overrides.
std::optional<std::string> showSystemParameterDialog(
    HWND owner, const std::string& system, const std::string& preset,
    const std::string& currentOverrides, const dsl::theme::Fonts& fonts);

}  // namespace dsl::compute
