#pragma once

#include "dynsyslab_viewer/Json.hpp"
#include "dynsyslab_viewer/Model.hpp"

#include <dynsyslab/core/Archive.hpp>

#include <filesystem>
#include <string>

namespace dsl {

[[nodiscard]] Payload buildLatticePayload(
    const std::filesystem::path& source, const std::string& view,
    const dynsyslab::core::NpzRecord& record, const json::Value& metadata);

}  // namespace dsl
