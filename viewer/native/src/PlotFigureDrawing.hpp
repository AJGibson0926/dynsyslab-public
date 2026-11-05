#pragma once
#include "dynsyslab_viewer/PlotCanvas.hpp"
namespace dsl {
void drawText(Gdiplus::Graphics&, const std::wstring&, const Gdiplus::RectF&, HFONT, COLORREF,
              Gdiplus::StringAlignment = Gdiplus::StringAlignmentNear);
PlotTransform renderTransform(const Figure&, const Gdiplus::RectF&, const std::optional<AxisLimits>&, const theme::Fonts&);
void figureCard(Gdiplus::Graphics&, const Figure&, const Gdiplus::RectF&,
    const theme::Fonts&, const PlotTransform&, const Plot3DCamera* = nullptr,
    Plot3DControl = Plot3DControl::none);
}
