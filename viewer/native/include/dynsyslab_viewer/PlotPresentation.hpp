#pragma once
#include "Model.hpp"
namespace dsl {
struct AxisLimits { double xMin = -1, xMax = 1, yMin = -1, yMax = 1; };
struct PhaseSpaceDomain {
    AxisLimits bounds;
    std::string coordinates;
    bool operator==(const PhaseSpaceDomain& other) const {
        return coordinates == other.coordinates && bounds.xMin == other.bounds.xMin
            && bounds.xMax == other.bounds.xMax && bounds.yMin == other.bounds.yMin
            && bounds.yMax == other.bounds.yMax;
    }
    bool operator!=(const PhaseSpaceDomain& other) const { return !(*this == other); }
};
bool plotReferenceVisible(const Figure& figure, const Figure::ReferencePoint& reference);
AxisLimits automaticAxisLimits(const Figure& figure);

bool canDisplayCoordinates(const Figure&,const std::string& coordinates);
void changeDisplayCoordinates(Figure&,const std::string& coordinates);
}
