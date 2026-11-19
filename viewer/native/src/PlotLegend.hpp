#pragma once

#include "dynsyslab_viewer/Model.hpp"
#include "dynsyslab_viewer/Theme.hpp"

#include <gdiplus.h>

namespace dsl {

// Draws the legend at `origin` (top-left corner), or at the default top-right corner.
void drawPlotLegend(Gdiplus::Graphics& graphics, const Figure& figure,
                    const Gdiplus::RectF& plot, const theme::Fonts& fonts,
                    const Gdiplus::PointF* origin = nullptr);
Gdiplus::SizeF plotLegendSize(Gdiplus::Graphics& graphics, const Figure& figure,
                              const Gdiplus::RectF& plot, const theme::Fonts& fonts);

}  // namespace dsl
