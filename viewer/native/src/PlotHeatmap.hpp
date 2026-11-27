#pragma once
#include "dynsyslab_viewer/PlotInteractionModel.hpp"
#include "dynsyslab_viewer/Theme.hpp"
namespace dsl {
void drawScalarHeatmap(Gdiplus::Graphics&, const Trace&, const PlotTransform&, const ValueScale& = {});
// The value range a heat map's colours span (fixed limits, else its finite values).
std::pair<double,double> heatmapColorRange(const Trace&, const ValueScale& = {});
}
