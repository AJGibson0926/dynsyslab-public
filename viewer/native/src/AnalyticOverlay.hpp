#pragma once
#include "DisplayCoordinates.hpp"
#include "dynsyslab_viewer/NativePayload.hpp"
#include <dynsyslab/core/Archive.hpp>
namespace dsl {
void appendAnalyticOverlay(Payload &payload, const dynsyslab::core::NpzRecord &record,
                           const DisplayCoordinateTransform &coordinates,
                           const std::vector<double> &time, const std::vector<std::size_t> &indices,
                           std::size_t samples, std::size_t dimension);
}
