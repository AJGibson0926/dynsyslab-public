#pragma once

#include "dynsyslab_viewer/Plot3D.hpp"
#include "dynsyslab_viewer/Theme.hpp"

#include <gdiplus.h>

namespace dsl {

void drawPlot3D(Gdiplus::Graphics& graphics, const Figure& figure,
                const Plot3DViewport& viewport, const Plot3DCamera& camera,
                Plot3DControl hoveredControl, const theme::Fonts& fonts);

}  // namespace dsl
