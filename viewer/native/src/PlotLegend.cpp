#include "PlotLegend.hpp"

#include "dynsyslab_viewer/MathText.hpp"
#include "MathTextDrawing.hpp"
#include "PlotSeries.hpp"

#include <algorithm>
#include <cmath>
#include <memory>

namespace dsl {
namespace {
struct LegendGeometry { std::vector<Trace> entries; std::size_t count = 0; float rowHeight = 0, sampleWidth = 0, sampleGap = 6.0f, width = 0; };
LegendGeometry legendGeometry(Gdiplus::Graphics& graphics, const Figure& figure,
                              const Gdiplus::RectF& plot, const theme::Fonts& fonts) {
    LegendGeometry g;
    for (const auto& trace : figure.traces) if (!trace.label.empty()) g.entries.push_back(trace);
    g.count = std::min<std::size_t>(g.entries.size(), 7);
    if (!g.count) return g;
    float textWidth = 0.0f, textHeight = 0.0f;
    for (std::size_t index = 0; index < g.count; ++index) {
        const auto size = measureMathText(graphics, wide(scientificMathLabel(g.entries[index].label)), fonts.caption);
        textWidth = std::max(textWidth, size.Width);
        textHeight = std::max(textHeight, size.Height);
    }
    g.rowHeight = std::max(24.0f, std::ceil(textHeight) + 8.0f);
    const bool hasTarget = std::any_of(g.entries.begin(),g.entries.begin()+g.count,
        [](const auto& trace) { return trace.role == "control_target" && trace.lineStyle == "dash"; });
    g.sampleWidth = std::max(hasTarget ? 50.0f : 20.0f, g.rowHeight * 0.88f);
    g.width = std::min(plot.Width * 0.62f, g.sampleWidth + g.sampleGap + textWidth + 3.0f);
    return g;
}
}  // namespace

Gdiplus::SizeF plotLegendSize(Gdiplus::Graphics& graphics, const Figure& figure,
                              const Gdiplus::RectF& plot, const theme::Fonts& fonts) {
    const auto g = legendGeometry(graphics, figure, plot, fonts);
    return {g.width, g.rowHeight * static_cast<float>(g.count)};
}

void drawPlotLegend(Gdiplus::Graphics& graphics, const Figure& figure,
                    const Gdiplus::RectF& plot, const theme::Fonts& fonts,
                    const Gdiplus::PointF* origin) {
    std::vector<Trace> entries;
    for (const auto& trace : figure.traces) if (!trace.label.empty()) entries.push_back(trace);
    const std::size_t count = std::min<std::size_t>(entries.size(), 7);
    if (!count) return;
    float textWidth = 0.0f, textHeight = 0.0f;
    for (std::size_t index = 0; index < count; ++index) {
        const auto size = measureMathText(graphics, wide(scientificMathLabel(entries[index].label)), fonts.caption);
        textWidth = std::max(textWidth, size.Width);
        textHeight = std::max(textHeight, size.Height);
    }
    const float rowHeight = std::max(24.0f, std::ceil(textHeight) + 8.0f);
    const bool hasTarget = std::any_of(entries.begin(),entries.begin()+count,
        [](const auto& trace) { return trace.role == "control_target" && trace.lineStyle == "dash"; });
    const float sampleWidth = std::max(hasTarget ? 50.0f : 20.0f, rowHeight * 0.88f);
    const float sampleGap = 6.0f;
    const float width = std::min(plot.Width * 0.62f,
        sampleWidth + sampleGap + textWidth + 3.0f);
    const Gdiplus::RectF area(origin ? origin->X : plot.GetRight() - width - 6.0f,
        origin ? origin->Y : plot.Y + 6.0f, width, rowHeight * static_cast<float>(count));
    for (std::size_t index = 0; index < count; ++index) {
        Gdiplus::Pen pen(plotSeriesColor(entries[index], index));
        stylePlotSeriesPen(pen,entries[index]);
        const float rowY = area.Y + static_cast<float>(index) * rowHeight;
        const float y = rowY + rowHeight * 0.5f;
        if (entries[index].encoding == "time") {
            for (int step = 0; step < 16; ++step) {
                pen.SetColor(applyPlotSeriesOpacity(trajectoryTimeColor(step / 15.),entries[index]));
                graphics.DrawLine(&pen, area.X + sampleWidth * step / 16, y,
                    area.X + sampleWidth * (step + 1) / 16, y);
            }
        } else graphics.DrawLine(&pen, area.X, y, area.X + sampleWidth, y);
        drawMathText(graphics, wide(scientificMathLabel(entries[index].label)),
            {area.X + sampleWidth + sampleGap, rowY,
             area.Width - sampleWidth - sampleGap, rowHeight},
            fonts.caption, theme::plotText);
    }
}

}  // namespace dsl
