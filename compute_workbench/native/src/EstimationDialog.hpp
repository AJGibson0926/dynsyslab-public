#pragma once

#include "dynsyslab_compute/Registry.hpp"
#include "dynsyslab_viewer/Theme.hpp"

#include <windows.h>

#include <filesystem>
#include <optional>
#include <string>

namespace dsl::compute {

// Configures observability analysis and KLQE from named reconstruction maps.
// The returned object is consumed unchanged by the canonical native workflow.
std::optional<std::string> showEstimationDialog(
    HWND owner, const std::filesystem::path& trajectorySource,
    const Registry& registry, const std::string& currentSettings,
    const dsl::theme::Fonts& fonts);

}  // namespace dsl::compute
