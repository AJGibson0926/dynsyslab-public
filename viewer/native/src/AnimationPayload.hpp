#pragma once
#include "dynsyslab_viewer/Model.hpp"

#include <dynsyslab/core/Archive.hpp>

#include <filesystem>

namespace dsl {

// Figures for a dynsyslab.animation/v1 manifest (AnimationPayload.cpp).
[[nodiscard]] Payload buildAnimationPayload(const std::filesystem::path& source, const dynsyslab::core::NpzRecord& record);

}  // namespace dsl
