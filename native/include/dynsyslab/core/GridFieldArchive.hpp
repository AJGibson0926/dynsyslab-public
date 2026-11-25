#pragma once
// dynsyslab.grid-field/v1: one or more GridFields on the same DomainGrid (frames, e.g. the Koopman
// action g o Phi^{k dt}), written by the domain_field job and read by Viewer.
//   arrays    axis_<j> (node values), anchor, field_values (frames x points x outputs, complex),
//             status (frames x points: 0 valid, 1 undefined, 2 conflict, 3 no preimage),
//             spread and counts (push-forwards), frame_time (frames), metadata_json
//   metadata  format, grid {domain, axes, coordinate_names}, outputs, output_names, complex,
//             frames, plus the producer's own keys (operation, function, map, source ...)
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/DomainGrid.hpp>
#include <dynsyslab/core/Json.hpp>

#include <vector>

namespace dynsyslab::core {

struct GridFieldSeries {
    std::vector<GridField> frames;
    std::vector<double> frameTime;  // one per frame (0 for a single field)
    json::Value metadata;
};

[[nodiscard]] NpzRecord gridFieldArchive(const GridFieldSeries& series);
[[nodiscard]] GridFieldSeries readGridFieldArchive(const NpzRecord& archive);
[[nodiscard]] bool isGridFieldArchive(const json::Value& metadata);

}  // namespace dynsyslab::core
