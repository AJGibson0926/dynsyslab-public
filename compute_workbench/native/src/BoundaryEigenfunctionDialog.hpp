#pragma once

#include "dynsyslab_compute/Registry.hpp"
#include "dynsyslab_viewer/Theme.hpp"

#include <windows.h>

#include <filesystem>
#include <optional>
#include <string>

namespace dsl::compute {

std::optional<std::string> showBoundaryEigenfunctionDialog(
    HWND owner,const std::filesystem::path& trajectorySource,
    const Registry& registry,const std::string& currentSettings,
    const dsl::theme::Fonts& fonts);

}  // namespace dsl::compute
