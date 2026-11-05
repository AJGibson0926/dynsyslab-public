#pragma once

#include "dynsyslab_viewer/Model.hpp"

#include <dynsyslab/core/Archive.hpp>

#include <filesystem>
#include <string>

namespace dsl {

// Builds the canonical spectrum plus one five-panel page per retained mode.
[[nodiscard]] Payload buildKromDashboardPayload(
    const std::filesystem::path& source,
    const dynsyslab::core::NpzRecord& archive,
    const std::string& metadataJson);

}  // namespace dsl
