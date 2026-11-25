#pragma once
// View "observable_phase": an observable's values along a trajectory (a time series of observable
// values, such as an eigenfunction evaluated by core::Observable::evaluate(TimeSeries)) drawn on that
// trajectory's own phase plane. The states at the evaluated times are coloured by the value: domain
// colouring for a complex observable (hue = argument, brightness = modulus), a colour scale for a
// real one. The whole source trajectory is drawn faintly underneath.
#include "dynsyslab_viewer/Json.hpp"
#include "dynsyslab_viewer/Model.hpp"

#include <dynsyslab/core/Archive.hpp>

#include <filesystem>

namespace dsl {

// True when the archive is a time series of observable values whose source trajectory is named.
[[nodiscard]] inline bool isObservableSeries(const json::Value& metadata) {
    return metadata.at("value_space").string() == "observable-values" && !metadata.at("source_archive").string().empty();
}

[[nodiscard]] Payload buildObservablePhasePayload(const std::filesystem::path& source,
                                                  const dynsyslab::core::NpzRecord& record, const json::Value& metadata);

}  // namespace dsl
