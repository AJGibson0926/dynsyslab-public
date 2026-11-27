#include "PlotLegend.hpp"

#include "dynsyslab_viewer/MathText.hpp"
#include "MathTextDrawing.hpp"
#include "PlotSeries.hpp"

#include "dynsyslab_viewer/PlotPresentation.hpp"
#include "dynsyslab_viewer/PlotScale.hpp"

#include <algorithm>
#include <array>
#include <optional>
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

namespace {
// One legend entry: its line sample, then its label; entries flow left to right and wrap.
struct TopLegendLayout {
    std::vector<Trace> entries;
    std::vector<float> widths;               // sample + gap + label
    std::vector<std::size_t> rowStarts;      // first entry of each row
    std::vector<float> rowWidths;
    float rowHeight = 0, sampleWidth = 0, sampleGap = 6.0f, entryGap = 26.0f;
};
TopLegendLayout topLegendLayout(Gdiplus::Graphics& graphics, const Figure& figure, float width,
                                const theme::Fonts& fonts) {
    TopLegendLayout g;
    for (const auto& trace : figure.traces) if (!trace.label.empty()) g.entries.push_back(trace);
    if (g.entries.empty()) return g;
    float textHeight = 0.0f;
    std::vector<float> text;
    for (const auto& entry : g.entries) {
        const auto size = measureMathText(graphics, wide(scientificMathLabel(entry.label)), fonts.caption);
        text.push_back(size.Width); textHeight = std::max(textHeight, size.Height);
    }
    g.rowHeight = std::max(24.0f, std::ceil(textHeight) + 8.0f);
    const bool hasTarget = std::any_of(g.entries.begin(), g.entries.end(),
        [](const auto& trace) { return trace.lineStyle == "dash"; });
    g.sampleWidth = std::max(hasTarget ? 50.0f : 20.0f, g.rowHeight * 0.88f);
    float row = 0.0f;
    for (std::size_t k = 0; k < g.entries.size(); ++k) {
        g.widths.push_back(g.sampleWidth + g.sampleGap + text[k] + 3.0f);
        const float next = row == 0.0f ? g.widths[k] : row + g.entryGap + g.widths[k];
        if (row == 0.0f || next > width) {
            if (row != 0.0f) g.rowWidths.push_back(row);
            g.rowStarts.push_back(k); row = g.widths[k];
        } else row = next;
    }
    g.rowWidths.push_back(row);
    return g;
}
}  // namespace

float plotTopLegendHeight(Gdiplus::Graphics& graphics, const Figure& figure, float width,
                          const theme::Fonts& fonts) {
    const auto g = topLegendLayout(graphics, figure, width, fonts);
    return g.rowHeight * static_cast<float>(g.rowStarts.size());
}

bool redundantPlotTitle(const Figure& figure) {
    return figure.title == "Phase space" || figure.title == "Control history"
        || (figure.kind == "time" && figure.title.find(" versus time") != std::string::npos);
}

std::optional<Gdiplus::RectF> clearLegendCorner(const Figure& figure, const Gdiplus::RectF& plot,
                                                const theme::Fonts& fonts) {
    if (figure.kind != "time" && figure.kind != "fourier") return std::nullopt;
    Gdiplus::Bitmap probe(1, 1);  // text measurement only
    Gdiplus::Graphics graphics(&probe);
    const auto size = plotLegendSize(graphics, figure, plot, fonts);
    if (!(size.Width > 0 && size.Height > 0)) return std::nullopt;
    const auto limits = automaticAxisLimits(figure);
    const auto valueScale = ordinateScale(figure);
    const double top = valueScale.forward(limits.yMax), span = top - valueScale.forward(limits.yMin);
    const double xSpan = limits.xMax - limits.xMin;
    const float gap = 6.0f, pad = 4.0f;  // the legend's offset from the frame, and a margin around it
    const std::array<Gdiplus::RectF, 4> corners{{
        {plot.GetRight() - size.Width - gap, plot.Y + gap, size.Width, size.Height},                       // top right
        {plot.X + gap, plot.Y + gap, size.Width, size.Height},                                             // top left
        {plot.GetRight() - size.Width - gap, plot.GetBottom() - size.Height - gap, size.Width, size.Height},  // bottom right
        {plot.X + gap, plot.GetBottom() - size.Height - gap, size.Width, size.Height}}};                   // bottom left
    for (auto corner : corners) {
        corner.Inflate(pad, pad);
        bool clear = true;
        for (const auto& trace : figure.traces) {
            if (trace.heatmap || !(xSpan > 0) || !(span != 0)) continue;
            const std::size_t count = std::min(trace.x.size(), trace.y.size());
            for (std::size_t i = 0; i < count && clear; ++i) {
                if (!std::isfinite(trace.x[i]) || !std::isfinite(trace.y[i])) continue;
                const float x = plot.X + static_cast<float>((trace.x[i] - limits.xMin) / xSpan) * plot.Width;
                const float y = plot.Y + static_cast<float>((top - valueScale.forward(trace.y[i])) / span) * plot.Height;
                clear = !corner.Contains(x, y);
            }
            if (!clear) break;
        }
        if (clear) { corner.Inflate(-pad, -pad); return corner; }
    }
    return std::nullopt;
}

// The band between the title (or the card's top edge) and the plot that holds the horizontal
// legend: always for "top"; for the automatic placement of a time or spectrum plot, only when the
// legend inside the plot (`plot`, laid out without a band) would cover data in every corner.
std::optional<Gdiplus::RectF> topLegendBand(const Figure& figure, const Gdiplus::RectF& card,
                                            const Gdiplus::RectF& plot, const theme::Fonts& fonts) {
    const bool automatic = figure.legendPlacement.empty() || figure.legendPlacement == "auto";
    if (!automatic && figure.legendPlacement != "top") return std::nullopt;
    if (automatic && (figure.kind != "time" && figure.kind != "fourier" || clearLegendCorner(figure, plot, fonts)))
        return std::nullopt;
    const float scale = scientificTextScale(fonts.body);
    const bool redundantTitle = redundantPlotTitle(figure);
    const float width = card.Width - 28.0f;
    Gdiplus::Bitmap probe(1, 1);
    Gdiplus::Graphics graphics(&probe);
    const float height = plotTopLegendHeight(graphics, figure, width, fonts);
    if (height <= 0.0f) return std::nullopt;
    return Gdiplus::RectF(card.X + 14.0f, card.Y + (redundantTitle ? 8.0f : 38.0f) * scale, width, height);
}

void drawPlotTopLegend(Gdiplus::Graphics& graphics, const Figure& figure,
                       const Gdiplus::RectF& band, const theme::Fonts& fonts) {
    const auto g = topLegendLayout(graphics, figure, band.Width, fonts);
    for (std::size_t row = 0; row < g.rowStarts.size(); ++row) {
        const std::size_t first = g.rowStarts[row];
        const std::size_t last = row + 1 < g.rowStarts.size() ? g.rowStarts[row + 1] : g.entries.size();
        float x = band.X + std::max(0.0f, (band.Width - g.rowWidths[row]) * 0.5f);
        const float rowY = band.Y + static_cast<float>(row) * g.rowHeight, y = rowY + g.rowHeight * 0.5f;
        for (std::size_t k = first; k < last; ++k) {
            const auto& entry = g.entries[k];
            Gdiplus::Pen pen(plotSeriesColor(entry, k));
            stylePlotSeriesPen(pen, entry);
            graphics.DrawLine(&pen, x, y, x + g.sampleWidth, y);
            drawMathText(graphics, wide(scientificMathLabel(entry.label)),
                {x + g.sampleWidth + g.sampleGap, rowY, g.widths[k] - g.sampleWidth - g.sampleGap + 2.0f, g.rowHeight},
                fonts.caption, theme::plotText);
            x += g.widths[k] + g.entryGap;
        }
    }
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
