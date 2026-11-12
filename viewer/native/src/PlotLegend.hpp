#pragma once

#include "dynsyslab_viewer/Model.hpp"
#include "dynsyslab_viewer/Theme.hpp"

#include <gdiplus.h>

namespace dsl {

void drawPlotLegend(Gdiplus::Graphics& graphics, const Figure& figure,
                    const Gdiplus::RectF& plot, const theme::Fonts& fonts);

}  // namespace dsl
