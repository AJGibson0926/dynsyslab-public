#pragma once

#include "dynsyslab_viewer/Theme.hpp"

#include <windows.h>

#include <optional>
#include <string>
#include <vector>

namespace dsl::compute {

// Edits a system-agnostic cartesian campaign. Each axis is a normal workflow
// JSON path, so the same UI covers parameter, forcing, target, time, solver,
// and initial-condition sweeps without embedding system-specific algorithms.
std::optional<std::string> showCampaignDialog(
    HWND owner, const std::string& current,
    const std::vector<std::string>& pathSuggestions,
    const dsl::theme::Fonts& fonts);

}  // namespace dsl::compute
