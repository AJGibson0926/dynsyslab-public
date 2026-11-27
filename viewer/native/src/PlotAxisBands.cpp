#include "PlotAxisBands.hpp"

#include "dynsyslab_viewer/MathText.hpp"
#include "MathTextDrawing.hpp"

#include <algorithm>
#include <cmath>
#include <cwchar>
#include <utility>

namespace dsl {

bool drawAxisBands(Gdiplus::Graphics& graphics, const Gdiplus::RectF& area, double yLow, double yHigh,
                   const Figure& figure, const theme::Fonts& fonts, const ValueScale& valueScale,
                   int yExponent, const std::wstring& leftTitle) {
    if (!figure.rightAxis || !figure.rightAxis->leftBand || !figure.rightAxis->rightBand) return false;
    const auto& axis = *figure.rightAxis;
    const float scale = scientificTextScale(fonts.body);
    const double top = valueScale.forward(yHigh), span = top - valueScale.forward(yLow);
    Gdiplus::Pen gridPen(theme::color(theme::plotGrid), 1.0f);
    // Four ticks across one band; returns the band's extent as fractions of the plot height.
    const auto side = [&](const std::pair<double, double>& band, bool right) {
        float upper = 1.0f, lower = 0.0f;
        for (int k = 0; k <= 3; ++k) {
            const double plotted = band.first + (band.second - band.first) * k / 3.0;
            const float fraction = static_cast<float>((top - valueScale.forward(plotted)) / span);
            if (!std::isfinite(fraction) || fraction < -0.01f || fraction > 1.01f) continue;
            upper = std::min(upper, fraction); lower = std::max(lower, fraction);
            const float y = area.Y + std::clamp(fraction, 0.0f, 1.0f) * area.Height;
            graphics.DrawLine(&gridPen, area.X, y, area.GetRight(), y);
            wchar_t label[40]{};
            swprintf_s(label, L"%.3g", right ? (plotted - axis.offset) / axis.scale : plotted / std::pow(10., yExponent));
            if (right) drawMathText(graphics, label, {area.GetRight() + 6*scale, y - 10*scale, 62*scale, 20*scale},
                                    fonts.body, theme::plotMuted, Gdiplus::StringAlignmentNear);
            else drawMathText(graphics, label, {area.X - 64*scale, y - 10*scale, 58*scale, 20*scale},
                              fonts.body, theme::plotMuted, Gdiplus::StringAlignmentFar);
        }
        return std::pair{std::min(upper, lower), std::max(upper, lower)};
    };
    const auto left = side(*axis.leftBand, false), right = side(*axis.rightBand, true);
    const auto band = [&](const std::pair<float, float>& extent, float x) {
        return Gdiplus::RectF(x, area.Y + extent.first * area.Height, 32*scale,
                              std::max(40*scale, (extent.second - extent.first) * area.Height));
    };
    drawVerticalText(graphics, leftTitle, band(left, area.X - 89*scale), fonts.body, theme::plotMuted);
    drawVerticalText(graphics, wide(scientificMathLabel(axis.label)), band(right, area.GetRight() + 68*scale),
                     fonts.body, theme::plotMuted);
    return true;
}

}  // namespace dsl
