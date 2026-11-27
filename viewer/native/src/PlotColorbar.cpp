#include "PlotColorbar.hpp"

#include "DomainColor.hpp"
#include "MathTextDrawing.hpp"
#include "PlotHeatmap.hpp"
#include "PlotSeries.hpp"
#include "dynsyslab_viewer/MathText.hpp"
#include "dynsyslab_viewer/PlotScale.hpp"

#include <algorithm>
#include <cstdio>

namespace dsl {
namespace {
constexpr double kPi = 3.14159265358979323846;
constexpr float kReserve = 74.0f;  // bar, gap and tick labels

bool domainEncoded(const Trace& trace) {
    return trace.encoding == "domain-color" || trace.encoding == "slice-domain-color";
}
bool scalarKind(const Figure& figure) {
    return figure.kind == "observable_heatmap" || figure.kind == "krom_heatmap" || figure.kind == "domain_field";
}
std::wstring number(double value) {
    wchar_t text[32]{};
    swprintf_s(text, L"%.3g", value);
    return text;
}
}  // namespace

std::optional<PlotColorbar> plotColorbar(const Figure& figure) {
    for (const auto& trace : figure.traces)
        if (domainEncoded(trace))
            return PlotColorbar{PlotColorbar::Kind::phase, -kPi, kPi,
                                figure.colorLabel.empty() ? std::string("\\arg") : figure.colorLabel};
    for (const auto& trace : figure.traces)
        if (trace.encoding == "time" && !trace.colorValues.empty()) {
            const auto [low, high] = trajectoryTimeRange(figure);
            return PlotColorbar{PlotColorbar::Kind::time, low, high, figure.colorLabel};
        }
    if (scalarKind(figure) && !figure.traces.empty()) {
        const auto [low, high] = heatmapColorRange(figure.traces.front(), resolveValueScale(figure));
        return PlotColorbar{PlotColorbar::Kind::scalar, low, high, figure.colorLabel};
    }
    return std::nullopt;
}

float plotColorbarReserve(const Figure& figure) { return plotColorbar(figure) ? kReserve : 0.0f; }

void drawPlotColorbar(Gdiplus::Graphics& graphics, const Figure& figure, const Gdiplus::RectF& plot,
                      const theme::Fonts& fonts) {
    const auto bar = plotColorbar(figure);
    if (!bar) return;
    const float scale = scientificTextScale(fonts.body);
    const float x = plot.GetRight() + 12 * scale, width = 12 * scale;
    const float top = plot.Y + 18 * scale, height = std::max(1.0f, plot.Height - 18 * scale);
    const auto valueScale = resolveValueScale(figure);
    const bool logScale = bar->kind == PlotColorbar::Kind::scalar && valueScale.kind != ValueScaleKind::linear;
    const double lower = logScale ? valueScale.forward(bar->low) : bar->low;
    const double upper = logScale ? valueScale.forward(bar->high) : bar->high;
    constexpr int slices = 128;
    for (int i = 0; i < slices; ++i) {  // bottom = low, top = high
        const double fraction = (i + 0.5) / slices;
        const auto color = bar->kind == PlotColorbar::Kind::phase ? domainColor(lower + (upper - lower) * fraction)
                         : bar->kind == PlotColorbar::Kind::time ? trajectoryTimeColor(fraction)
                                                                 : theme::heatColor(fraction);
        Gdiplus::SolidBrush brush(color);
        graphics.FillRectangle(&brush, x, top + height * (1.0f - float(i + 1) / slices), width, height / slices + 0.5f);
    }
    Gdiplus::Pen outline(theme::color(theme::plotBorder), 1.0f);
    graphics.DrawRectangle(&outline, x, top, width, height);
    const auto tick = [&](double fraction, const std::wstring& text) {
        const float y = top + height * float(1.0 - fraction);
        graphics.DrawLine(&outline, x + width, y, x + width + 3 * scale, y);
        drawMathText(graphics, text, {x + width + 5 * scale, y - 10 * scale, kReserve * scale, 20 * scale},
                     fonts.caption, theme::plotMuted);
    };
    if (bar->kind == PlotColorbar::Kind::phase) {
        tick(0.0, wide(scientificMathLabel("-\\pi"))); tick(0.5, L"0"); tick(1.0, wide(scientificMathLabel("\\pi")));
    } else {
        const auto valueAt = [&](double fraction) {
            const double v = lower + (upper - lower) * fraction;
            return logScale ? valueScale.inverse(v) : v;
        };
        for (const double fraction : {0.0, 0.5, 1.0}) tick(fraction, number(valueAt(fraction)));
    }
    if (!bar->label.empty())  // what the colours mean, above the bar
        drawMathText(graphics, wide(scientificMathLabel(bar->label)),
                     {x - 4 * scale, plot.Y - 4 * scale, (kReserve + 8) * scale, 20 * scale}, fonts.caption,
                     theme::plotMuted);
}

}  // namespace dsl
