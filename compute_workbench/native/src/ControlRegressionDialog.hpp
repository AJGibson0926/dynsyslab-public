#pragma once

#include "dynsyslab_compute/Registry.hpp"
#include "dynsyslab_viewer/Theme.hpp"

#include <windows.h>

#include <filesystem>
#include <optional>
#include <string>

namespace dsl::compute {

// Configures constant and bilinear B_Psi regression against a controlled
// trajectory archive. The generated error-map grid is native Viewer data.
std::optional<std::string> showControlRegressionDialog(
    HWND owner, const std::filesystem::path& controlledSource,
    const Registry& registry, const std::string& currentSettings,
    const dsl::theme::Fonts& fonts);

}  // namespace dsl::compute
