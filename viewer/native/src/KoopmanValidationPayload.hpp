#pragma once

#include "dynsyslab_viewer/Model.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Json.hpp>

#include <filesystem>
#include <string>

namespace dsl {

Payload buildKoopmanValidationPayload(
    const std::filesystem::path& source,const std::string& view,
    const dynsyslab::core::NpzRecord& record,
    const dynsyslab::core::json::Value& metadata);

}  // namespace dsl
