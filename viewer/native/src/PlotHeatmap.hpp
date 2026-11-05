#pragma once
#include "dynsyslab_viewer/PlotInteractionModel.hpp"
#include "dynsyslab_viewer/Theme.hpp"
namespace dsl {
void drawScalarHeatmap(Gdiplus::Graphics&, const Trace&, const PlotTransform&, const ValueScale& = {});
void drawScalarColorbar(Gdiplus::Graphics&, const Figure&, const Gdiplus::RectF&, const theme::Fonts&);
}
