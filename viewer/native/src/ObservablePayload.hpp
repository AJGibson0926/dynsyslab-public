#pragma once

#include "dynsyslab_viewer/Json.hpp"
#include "dynsyslab_viewer/PlotPresentation.hpp"

#include <dynsyslab/core/Archive.hpp>

#include <filesystem>
#include <string>

namespace dsl {

[[nodiscard]] Payload buildObservablePayload(
    const std::filesystem::path& source, const std::string& view,
    const dynsyslab::core::NpzRecord& record, const json::Value& metadata,
    const std::optional<PhaseSpaceDomain>& phaseDomain = std::nullopt);

void appendObservablePayload(Payload& target, Payload source, const std::string& label,
    const std::string& color = {}, bool dashed = false);

}  // namespace dsl
