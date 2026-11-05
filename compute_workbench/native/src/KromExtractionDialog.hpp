#pragma once

#include "dynsyslab_viewer/Theme.hpp"

#include <windows.h>

#include <filesystem>
#include <optional>
#include <string>

namespace dsl::compute {

// Selects eigenfunctions from a self-contained KROM and returns settings for
// an independent extracted archive. Displayed indices are the source KROM's
// zero-based local indices; inherited ancestral IDs remain archive evidence.
std::optional<std::string> showKromExtractionDialog(
    HWND owner, const std::filesystem::path& source,
    const std::string& currentSettings, const dsl::theme::Fonts& fonts);

}  // namespace dsl::compute
