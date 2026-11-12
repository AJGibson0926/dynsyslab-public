#pragma once
#include "dynsyslab_viewer/NativePayload.hpp"
namespace dsl {
Payload buildSindyDynamicsPayload(const std::filesystem::path& source,const std::string& view,
    const dynsyslab::core::NpzRecord& record,const json::Value& metadata);
}
