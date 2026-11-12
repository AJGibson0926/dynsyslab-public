#pragma once

#include "dynsyslab_viewer/Model.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Json.hpp>

#include <filesystem>
#include <string>

namespace dsl {

Payload buildFourierPayload(const std::filesystem::path& source,
                            const dynsyslab::core::NpzRecord& archive,
                            const dynsyslab::core::json::Value& metadata);

Payload buildControlFourierPayload(
    const std::filesystem::path& source,
    const dynsyslab::core::NpzRecord& archive,
    const dynsyslab::core::json::Value& metadata);

}  // namespace dsl
