#pragma once

#include "dynsyslab_viewer/Model.hpp"
#include "dynsyslab_viewer/PlotInteractionModel.hpp"
#include "dynsyslab_viewer/Theme.hpp"

#include <gdiplus.h>

namespace dsl {

void drawCurveProbe(Gdiplus::Graphics& graphics, const Figure& figure,
                    const PlotTransform& transform, const CurveProbe& probe,
                    bool pinned, const theme::Fonts& fonts);

}  // namespace dsl
