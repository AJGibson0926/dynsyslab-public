#include "PlotProbeOverlay.hpp"
#include "PlotSeries.hpp"

#include <algorithm>
#include <iomanip>
#include <sstream>

namespace dsl {
namespace {

std::wstring coordinates(const CurveProbe& probe) {
    std::wostringstream text;
    text << std::setprecision(6) << L"(" << probe.x << L", " << probe.y;
    if (probe.z) text << L", " << *probe.z;
    text << L")";
    return text.str();
}

void tooltip(Gdiplus::Graphics& graphics, const PlotTransform& transform,
             const CurveProbe& probe, const theme::Fonts& fonts) {
    constexpr float height = 28.0f, gap = 10.0f;
    const float desiredWidth = probe.z ? 240.0f : 178.0f;
    const float width = std::max(1.0f,
        std::min(desiredWidth, static_cast<float>(transform.width)));
    float left = static_cast<float>(probe.pixelX) + gap;
    float top = static_cast<float>(probe.pixelY) - height - gap;
    const float right = static_cast<float>(transform.left + transform.width);
    const float bottom = static_cast<float>(transform.top + transform.height);
    if (left + width > right) left = static_cast<float>(probe.pixelX) - width - gap;
    if (top < transform.top) top = static_cast<float>(probe.pixelY) + gap;
    left = std::clamp(left, static_cast<float>(transform.left),
                      std::max(static_cast<float>(transform.left), right - width));
    top = std::clamp(top, static_cast<float>(transform.top),
                     std::max(static_cast<float>(transform.top), bottom - height));
    Gdiplus::SolidBrush fill(theme::color(theme::panelRaised, 242));
    Gdiplus::Pen outline(theme::color(theme::accent, 220), 1.0f);
    const Gdiplus::RectF area(left, top, width, height);
    graphics.FillRectangle(&fill, area); graphics.DrawRectangle(&outline, area);
    HDC dc = graphics.GetHDC(); Gdiplus::Font font(dc, fonts.caption);
    graphics.ReleaseHDC(dc);
    Gdiplus::SolidBrush textBrush(theme::color(theme::text));
    Gdiplus::StringFormat format;
    format.SetAlignment(Gdiplus::StringAlignmentCenter);
    format.SetLineAlignment(Gdiplus::StringAlignmentCenter);
    const auto label = coordinates(probe);
    graphics.DrawString(label.c_str(), -1, &font, area, &format, &textBrush);
}

}  // namespace

void drawCurveProbe(Gdiplus::Graphics& graphics, const Figure& figure,
                    const PlotTransform& transform, const CurveProbe& probe,
                    bool pinned, const theme::Fonts& fonts) {
    const auto color = probe.traceIndex < figure.traces.size()
        ? plotSeriesColor(figure.traces[probe.traceIndex], probe.traceIndex)
        : theme::color(theme::accent);
    const float x = static_cast<float>(probe.pixelX);
    const float y = static_cast<float>(probe.pixelY);
    if (pinned) {
        Gdiplus::SolidBrush halo(theme::color(theme::plotBackground, 235));
        graphics.FillEllipse(&halo, x - 6.0f, y - 6.0f, 12.0f, 12.0f);
    }
    Gdiplus::SolidBrush fill(color);
    Gdiplus::Pen ring(theme::color(pinned ? theme::plotText : theme::plotBackground), 1.4f);
    const float radius = pinned ? 4.2f : 3.5f;
    graphics.FillEllipse(&fill, x - radius, y - radius, radius * 2, radius * 2);
    graphics.DrawEllipse(&ring, x - radius, y - radius, radius * 2, radius * 2);
    if (pinned) tooltip(graphics, transform, probe, fonts);
}

}  // namespace dsl
