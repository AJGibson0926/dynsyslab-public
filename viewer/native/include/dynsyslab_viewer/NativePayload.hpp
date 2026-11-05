#pragma once

#include "dynsyslab_viewer/Model.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Json.hpp>

#include <filesystem>
#include <string>

namespace dsl {

// Converts a persisted scientific archive directly into lightweight Viewer traces.
// This service is deliberately independent of the windowing layer.
[[nodiscard]] Payload buildNativePayload(const std::filesystem::path& source,
                                         const std::string& view);

// Multi-panel clients may decode one archive once and request several
// independent views without repeating ZIP inflation.
[[nodiscard]] Payload buildNativePayload(
    const std::filesystem::path& source, const std::string& view,
    const dynsyslab::core::NpzRecord& record,
    const dynsyslab::core::json::Value& metadata);

}  // namespace dsl
