#pragma once

#include "dynsyslab_compute/Registry.hpp"
#include "dynsyslab_viewer/Theme.hpp"

#include <dynsyslab/core/LegacyCsv.hpp>

#include <windows.h>

#include <filesystem>
#include <optional>

namespace dsl::compute {

struct LegacyImportRequest {
    std::filesystem::path source;
    dynsyslab::core::LegacyCsvImportOptions options;
};

// Builds an explicit, reviewable mapping from a tagged legacy MATLAB CSV to
// canonical native trajectory archives. No legacy dynamics name is silently
// treated as an executable native model.
std::optional<LegacyImportRequest> showLegacyImportDialog(
    HWND owner, const Registry& registry, const dsl::theme::Fonts& fonts);

}  // namespace dsl::compute
