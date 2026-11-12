#pragma once

#include "dynsyslab_viewer/Theme.hpp"

#include <gdiplus.h>
#include <windows.h>

namespace dsl {

RECT panelRectangle(const Gdiplus::RectF& area);
void drawPanelChrome(Gdiplus::Graphics& graphics, const Gdiplus::RectF& card,
                     bool highlighted, const theme::Fonts& fonts);

}  // namespace dsl
