#pragma once
#include "DisplayCoordinates.hpp"
#include "dynsyslab_viewer/Model.hpp"
namespace dsl {
void decorateTrajectory(Payload&,const json::Value&,const DisplayCoordinateTransform&,
    const std::vector<double>& time);
}
