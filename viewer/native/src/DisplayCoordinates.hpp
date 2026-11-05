#pragma once

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab_viewer/Json.hpp>

#include <cstddef>
#include <string>
#include <vector>

namespace dsl {

struct DisplayCoordinateTransform {
    std::string coordinates;
    std::string timeLabel{"time"};
    std::string suffix;
    std::vector<std::string> stateLabels;
    double timeMultiplier{1.0};
    std::vector<double> stateMultipliers;
    std::string warning;
};

DisplayCoordinateTransform displayCoordinateTransform(
    const json::Value& metadata, std::size_t stateDimension,
    const std::string& requested = {});
std::vector<double> transformDisplayTime(
    std::vector<double> values, const DisplayCoordinateTransform& transform);
dynsyslab::core::NpyArray transformDisplayStates(
    const dynsyslab::core::NpyArray& values,
    const DisplayCoordinateTransform& transform);

}  // namespace dsl
