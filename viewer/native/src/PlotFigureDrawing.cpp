#include "dynsyslab_viewer/PlotCanvas.hpp"
#include "dynsyslab_viewer/PanelWorkspace.hpp"
#include "dynsyslab_viewer/PlotInteractionModel.hpp"
#include "PlotPanelChrome.hpp"
#include "PlotLegend.hpp"
#include "PlotSeries.hpp"
#include "PlotHeatmap.hpp"
#include "MathTextDrawing.hpp"
#include "dynsyslab_viewer/MathText.hpp"
#include "dynsyslab_viewer/PlotPresentation.hpp"
#include "PlotProbeOverlay.hpp"
#include "Plot3DDrawing.hpp"

#include <gdiplus.h>

#include <algorithm>
#include <cmath>
#include <limits>

#include "PlotFigureDrawing.hpp"
namespace dsl {
struct Range {
    double low = std::numeric_limits<double>::infinity();
    double high = -std::numeric_limits<double>::infinity();
    void include(double value) {
        if (std::isfinite(value)) { low = std::min(low, value); high = std::max(high, value); }
    }
    void normalize() {
        if (!std::isfinite(low) || !std::isfinite(high)) { low = -1; high = 1; }
        if (std::abs(high - low) < 1e-14) { low -= 0.5; high += 0.5; }
        const double margin = (high - low) * 0.045;
        low -= margin; high += margin;
    }
};
void drawText(Gdiplus::Graphics& graphics, const std::wstring& value,
              const Gdiplus::RectF& area, HFONT font, COLORREF color,
              Gdiplus::StringAlignment align) {
    HDC dc = graphics.GetHDC();
    Gdiplus::Font display(dc, font);
    graphics.ReleaseHDC(dc);
    Gdiplus::SolidBrush brush(theme::color(color));
    Gdiplus::StringFormat format;
    format.SetAlignment(align);
    format.SetLineAlignment(Gdiplus::StringAlignmentCenter);
    format.SetTrimming(Gdiplus::StringTrimmingEllipsisCharacter);
    graphics.DrawString(value.c_str(), -1, &display, area, &format, &brush);
}
void drawVerticalText(Gdiplus::Graphics& graphics, const std::wstring& value,
                      const Gdiplus::RectF& area, HFONT font, COLORREF color) {
    const Gdiplus::GraphicsState state = graphics.Save();
    graphics.TranslateTransform(area.X, area.GetBottom());
    graphics.RotateTransform(-90.0f);
    drawMathText(graphics, value, {0.0f, 0.0f, area.Height, area.Width}, font, color,
             Gdiplus::StringAlignmentCenter);
    graphics.Restore(state);
}
void rounded(Gdiplus::Graphics& graphics, const Gdiplus::RectF& area,
             COLORREF fillColor, COLORREF outlineColor, float radius = 10.0f) {
    Gdiplus::GraphicsPath path;
    path.AddArc(area.X, area.Y, radius, radius, 180, 90);
    path.AddArc(area.GetRight() - radius, area.Y, radius, radius, 270, 90);
    path.AddArc(area.GetRight() - radius, area.GetBottom() - radius, radius, radius, 0, 90);
    path.AddArc(area.X, area.GetBottom() - radius, radius, radius, 90, 90);
    path.CloseFigure();
    Gdiplus::SolidBrush fillBrush(theme::color(fillColor));
    Gdiplus::Pen outline(theme::color(outlineColor), 1.0f);
    graphics.FillPath(&fillBrush, &path);
    graphics.DrawPath(&outline, &path);
}
COLORREF colorRef(const std::string& encoded, std::size_t fallback) {
    const auto value = theme::plotColor(encoded, fallback);
    return RGB(value.GetR(), value.GetG(), value.GetB());
}
Gdiplus::RectF squareInside(const Gdiplus::RectF& value) {
    const float side = std::min(value.Width, value.Height);
    return {value.X + (value.Width - side) / 2.0f,
            value.Y + (value.Height - side) / 2.0f, side, side};
}
Gdiplus::RectF plotArea(const Figure& figure, const Gdiplus::RectF& card, const theme::Fonts& fonts) {
    const float scale = scientificTextScale(fonts.body);
    const bool redundantTitle = figure.title == "Phase space" ||
        figure.title == "Control history" ||
        (figure.kind == "time" && figure.title.find(" versus time") != std::string::npos);
    const float top = (figure.kind == "observable_heatmap" || figure.kind == "krom_heatmap" || figure.kind == "krom_time" && !figure.colorLabel.empty()) ? 70.0f : !figure.colorLabel.empty() ? 36.0f : (redundantTitle ? 10.0f : 44.0f);
    Gdiplus::RectF plot(card.X + 90*scale, card.Y + top*scale,
        std::max(1.0f, card.Width - 130*scale), std::max(1.0f, card.Height - (top+62)*scale));
    return figure.square && figure.kind != "phase" ? squareInside(plot) : plot;
}

void ranges(const Figure& figure, Range& xr, Range& yr) {
    const auto limits = automaticAxisLimits(figure);
    xr = {limits.xMin, limits.xMax}; yr = {limits.yMin, limits.yMax};
}

PlotTransform renderTransform(const Figure& figure, const Gdiplus::RectF& card,
                              const std::optional<AxisLimits>& limits, const theme::Fonts& fonts) {
    const auto plot = plotArea(figure, card, fonts);
    Range xr, yr; ranges(figure, xr, yr);
    if (limits) {
        xr = {limits->xMin, limits->xMax};
        yr = {limits->yMin, limits->yMax};
    }
    return {plot.X, plot.Y, plot.Width, plot.Height,
            xr.low, xr.high, yr.low, yr.high, ordinateScale(figure)};
}


void grid(Gdiplus::Graphics& graphics, const Gdiplus::RectF& area,
          const Range& xr, const Range& yr, const Figure& figure,
          const theme::Fonts& fonts) {
    const float scale = scientificTextScale(fonts.body);
    auto exponent = [](const Range& range) {
        const auto magnitude = std::max(std::abs(range.low), std::abs(range.high));
        return magnitude > 0 && (magnitude >= 10000 || magnitude < .001)
            ? static_cast<int>(std::floor(std::log10(magnitude))) : 0;
    };
    const auto valueScale = ordinateScale(figure);
    const int xExponent = exponent(xr), yExponent = valueScale.kind == ValueScaleKind::linear ? exponent(yr) : 0;
    auto axisLabel = [](const std::string& name, int power) {
        return wide(scientificMathLabel(name) + (power
            ? "\\; \\times 10^{" + std::to_string(power) + "}" : ""));
    };
    Gdiplus::Pen gridPen(theme::color(theme::plotGrid), 1.0f);
    for (int index = 0; index <= 4; ++index) {
        const float fraction = static_cast<float>(index) / 4.0f;
        graphics.DrawLine(&gridPen, area.X + fraction * area.Width, area.Y,
                          area.X + fraction * area.Width, area.GetBottom());
        graphics.DrawLine(&gridPen, area.X, area.Y + fraction * area.Height,
                          area.GetRight(), area.Y + fraction * area.Height);
        const double x = xr.low + fraction * (xr.high - xr.low);
        const double y = valueScale.inverse(valueScale.forward(yr.high) - fraction *
            (valueScale.forward(yr.high) - valueScale.forward(yr.low)));
        wchar_t label[40]{};
        swprintf_s(label, L"%.3g", x / std::pow(10., xExponent));
        drawMathText(graphics, label, {area.X + fraction * area.Width - 38*scale,
            area.GetBottom() + 2*scale, 76*scale, 20*scale}, fonts.body, theme::plotMuted,
            Gdiplus::StringAlignmentCenter);
        swprintf_s(label, L"%.3g", y / std::pow(10., yExponent));
        drawMathText(graphics, label, {area.X - 64*scale, area.Y + fraction * area.Height - 10*scale,
            58*scale, 20*scale}, fonts.body, theme::plotMuted, Gdiplus::StringAlignmentFar);
    }
    drawMathText(graphics, axisLabel(figure.xLabel, xExponent), {area.X, area.GetBottom() + 23*scale,
        area.Width, 34*scale}, fonts.body, theme::plotMuted, Gdiplus::StringAlignmentCenter);
    drawVerticalText(graphics, axisLabel(figure.yLabel + (valueScale.kind == ValueScaleKind::linear ? "" : " [" + valueScale.name() + "]"), yExponent), {area.X - 89*scale, area.Y,
        32*scale, area.Height}, fonts.body, theme::plotMuted);
}

void heatmap(Gdiplus::Graphics& graphics, const Trace& trace,
             const PlotTransform& view, const AxisLimits& full) {
    if (trace.z.empty()) return;
    int rows = 0, columns = 0;
    if (trace.zShape.size() >= 2) {
        rows = trace.zShape[trace.zShape.size() - 2];
        columns = trace.zShape.back();
    }
    if (rows <= 0 || columns <= 0 || static_cast<std::size_t>(rows * columns) > trace.z.size()) return;
    Range values;
    for (double value : trace.z) values.include(value);
    values.normalize();
    // Legacy matrices fill their automatic axes. Project those same cell
    // boundaries through the current view so the pixels follow zoom and pan.
    const auto xEdge = [&](int column) {
        const double x = full.xMin + (full.xMax - full.xMin) * (double(column) / columns);
        return static_cast<float>(view.left + view.width * std::clamp(
            (x - view.xMin) / (view.xMax - view.xMin), 0.0, 1.0));
    };
    const auto yEdge = [&](int row) {
        const double y = full.yMin + (full.yMax - full.yMin) * (double(row) / rows);
        return static_cast<float>(std::clamp(plotPixelY(view, y), view.top, view.top + view.height));
    };
    for (int row = 0; row < rows; ++row) {
        const float top = yEdge(row + 1), bottom = yEdge(row);
        if (!(bottom > top)) continue;
        for (int column = 0; column < columns; ++column) {
            const float left = xEdge(column), right = xEdge(column + 1);
            if (!(right > left)) continue;
            const double value = trace.z[static_cast<std::size_t>(row * columns + column)];
            const double normalized = (value - values.low) / (values.high - values.low);
            Gdiplus::SolidBrush brush(theme::heatColor(normalized));
            graphics.FillRectangle(&brush, left, top, right - left + 0.8f, bottom - top + 0.8f);
        }
    }
}

void inset(Gdiplus::Graphics& graphics, const Figure& figure,
           const Gdiplus::RectF& plot, const theme::Fonts& fonts) {
    const auto resolved = figure.inset ? figure.inset : ordinateScale(figure).kind == ValueScaleKind::linear ? automaticDetailInset(figure) : std::nullopt;
    if (!resolved) return;
    const Gdiplus::RectF area(plot.GetRight() - plot.Width * 0.39f - 8,
        plot.GetBottom() - plot.Height * 0.37f - 8, plot.Width * 0.39f, plot.Height * 0.37f);
    Gdiplus::SolidBrush background(theme::color(theme::plotBackground, 238));
    Gdiplus::Pen outline(theme::color(theme::accent, 190), 1.0f);
    graphics.FillRectangle(&background, area);
    graphics.DrawRectangle(&outline, area);
    const Gdiplus::RectF inner(area.X + 7, area.Y + 18, area.Width - 14, area.Height - 25);
    Range xr{resolved->xMin, resolved->xMax};
    Range yr{resolved->yMin, resolved->yMax};
    drawText(graphics, resolved->automatic ? L"automatic detail" : L"detail",
             {area.X + 7, area.Y, area.Width - 14, 18},
             fonts.caption, theme::accentStrong);
    for (std::size_t index = 0; index < figure.traces.size(); ++index) {
        if (!resolved->traceIndices.empty() && std::find(resolved->traceIndices.begin(),
                resolved->traceIndices.end(), index) == resolved->traceIndices.end()) continue;
        drawPlotSeries(graphics, figure.traces[index], index,
            {inner.X, inner.Y, inner.Width, inner.Height, xr.low, xr.high, yr.low, yr.high, ordinateScale(figure)});
    }
}

void figureCard(Gdiplus::Graphics& graphics, const Figure& figure,
                const Gdiplus::RectF& card, const theme::Fonts& fonts,
                const PlotTransform& transform, const Plot3DCamera* camera3D,
                Plot3DControl hovered3DControl) {
    rounded(graphics, card, theme::plotBackground, theme::plotBorder, 13);
    const bool redundantTitle = figure.title == "Phase space" ||
        figure.title == "Control history" ||
        (figure.kind == "time" && figure.title.find(" versus time") != std::string::npos);
    const bool showTitle = !redundantTitle;
    if (showTitle)
        drawMathText(graphics, wide(scientificMathLabel(figure.title)), {card.X + 14, card.Y + 5,
            card.Width - 28, 30}, fonts.heading,
            figure.titleColor.empty() ? theme::plotText : colorRef(figure.titleColor,0));
    if (!figure.traces.empty() && figure.traces.front().encoding == "domain-color")
        drawText(graphics,figure.kind == "krom_time" ? L"Hue: phase; height: magnitude" : L"Hue: phase; brightness: relative magnitude",{card.X+14,card.Y+31,card.Width-28,20},fonts.caption,theme::plotMuted);
    const Gdiplus::RectF plot(static_cast<float>(transform.left),
        static_cast<float>(transform.top), static_cast<float>(transform.width),
        static_cast<float>(transform.height));
    if (figure.kind == "three_d") {
        const Plot3DCamera fallback;
        drawPlot3D(graphics, figure,
            {transform.left, transform.top, transform.width, transform.height},
            camera3D ? *camera3D : fallback, hovered3DControl, fonts);
        drawTimeColorbar(graphics, figure, plot, fonts);
        return;
    }
    Range xr{transform.xMin, transform.xMax};
    Range yr{transform.yMin, transform.yMax};
    grid(graphics, plot, xr, yr, figure, fonts);
    // On phase portraits the time-colored solution stays visible over a closely tracking target.
    for (int pass = 0; pass < 3; ++pass)
        for (std::size_t index = 0; index < figure.traces.size(); ++index) {
            const Trace& trace = figure.traces[index];
            const bool coloredPhase = figure.kind == "phase" && trace.encoding == "time";
            const int layer = trace.role.rfind("phase_",0)==0 ? 0 : coloredPhase ? 2 : 1;
            if (layer != pass) continue;
            if (trace.heatmap && (figure.kind == "observable_heatmap" || figure.kind == "krom_heatmap")) drawScalarHeatmap(graphics,trace,transform,resolveValueScale(figure));
            else if (trace.heatmap) heatmap(graphics, trace, transform, automaticAxisLimits(figure));
            else drawPlotSeries(graphics, trace, index, transform, trajectoryTimeRange(figure));
        }
    drawPlotReferences(graphics, figure, transform, fonts);
    drawTimeColorbar(graphics, figure, plot, fonts);
    drawScalarColorbar(graphics, figure, plot, fonts);
    if (figure.kind != "krom_spectrum" && figure.kind != "stem"
        && !(figure.kind == "krom_time" && figure.traces.size() == 1)
        && !figure.traces.empty() && !figure.traces.front().heatmap)
        drawPlotLegend(graphics, figure, plot, fonts);
    if (figure.kind == "krom_spectrum") {
        Figure selected = figure; selected.traces.clear();
        for (const auto& t : figure.traces) if (t.selected) selected.traces.push_back(t);
        drawPlotLegend(graphics,selected,plot,fonts);
    }
    inset(graphics, figure, plot, fonts);
}

}
