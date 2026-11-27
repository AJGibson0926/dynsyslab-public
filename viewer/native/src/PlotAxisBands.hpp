#pragma once
#include "dynsyslab_viewer/Model.hpp"
#include "dynsyslab_viewer/PlotScale.hpp"
#include "dynsyslab_viewer/Theme.hpp"

#include <gdiplus.h>
#include <string>

namespace dsl {
void drawVerticalText(Gdiplus::Graphics& graphics, const std::wstring& value,
                      const Gdiplus::RectF& area, HFONT font, COLORREF color);
// Stacked ordinates (Figure::RightAxis leftBand/rightBand): the left curves' band with its ticks,
// grid lines and title on the left, the right curves' band with its own on the right. Returns
// false, drawing nothing, for a figure without bands. Left tick values are divided by 10^yExponent.
bool drawAxisBands(Gdiplus::Graphics& graphics, const Gdiplus::RectF& area, double yLow, double yHigh,
                   const Figure& figure, const theme::Fonts& fonts, const ValueScale& valueScale,
                   int yExponent, const std::wstring& leftTitle);
}  // namespace dsl
