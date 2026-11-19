#pragma once
#include "dynsyslab_viewer/Model.hpp"
namespace dsl {
enum class ValueScaleKind { linear, logarithmic, symmetric_log };
struct ValueScale {
    ValueScaleKind kind = ValueScaleKind::linear;
    double linearThreshold = 1;
    std::string reason;
    // Symmetric log only: plotted height of one decade beyond the linear core,
    // relative to the core half-height (1 = standard symlog).
    double decadeWeight = 1;
    double forward(double value) const;
    double inverse(double value) const;
    std::string name() const;
};
bool supportsValueScale(const Figure& figure);
bool scalesColorValues(const Figure& figure);
ValueScale resolveValueScale(const Figure& figure);
ValueScale ordinateScale(const Figure& figure);
void setValueScale(Payload& payload, const std::string& requested);
}
