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
#include "dynsyslab_viewer/PlotScale.hpp"
#include "PlotProbeOverlay.hpp"
#include "Plot3DDrawing.hpp"

#include <gdiplus.h>

#include <algorithm>
#include <array>
#include <vector>
#include <optional>
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
    const auto xTick = [&](double fraction, double x) {
        graphics.DrawLine(&gridPen, area.X + static_cast<float>(fraction) * area.Width, area.Y,
                          area.X + static_cast<float>(fraction) * area.Width, area.GetBottom());
        wchar_t label[40]{};
        swprintf_s(label, L"%.3g", x / std::pow(10., xExponent));
        drawMathText(graphics, label, {area.X + static_cast<float>(fraction) * area.Width - 38*scale,
            area.GetBottom() + 2*scale, 76*scale, 20*scale}, fonts.body, theme::plotMuted,
            Gdiplus::StringAlignmentCenter);
    };
    const double first = std::ceil(xr.low), last = std::floor(xr.high);
    const bool integers = figure.integerXTicks && xr.high > xr.low && last >= first && last - first <= 1e4;
    if (integers) {  // one tick per integer, thinned to at most ~12 labels
        const double step = std::max(1.0, std::ceil((last - first + 1) / 12.0));
        for (double x = first; x <= last; x += step) xTick((x - xr.low) / (xr.high - xr.low), x);
    }
    for (int index = 0; index <= 4; ++index) {
        const float fraction = static_cast<float>(index) / 4.0f;
        if (!integers) xTick(fraction, xr.low + fraction * (xr.high - xr.low));
        graphics.DrawLine(&gridPen, area.X, area.Y + fraction * area.Height,
                          area.GetRight(), area.Y + fraction * area.Height);
        const double y = valueScale.inverse(valueScale.forward(yr.high) - fraction *
            (valueScale.forward(yr.high) - valueScale.forward(yr.low)));
        wchar_t label[40]{};
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

// Overlay placement. The legend sits top right and the detail inset bottom right, and
// either moves only for a real reason: it would hide something that matters in the main
// panel, namely each trace's visible extremes or the episode the inset magnifies.
// Ordinary repetitive data under an overlay is not a reason to move it.
struct OverlayLayout { std::optional<Gdiplus::RectF> legend, inset; double cost = 0; };
OverlayLayout overlayLayout(const Figure& figure, const Gdiplus::RectF& plot, const PlotTransform& view,
                            std::optional<Gdiplus::SizeF> legendSize, std::optional<Gdiplus::SizeF> insetSize,
                            const std::optional<DetailInset>& detail) {
    auto corners = [&](float width, float height, float gap) {
        return std::array<Gdiplus::RectF, 4>{{
            {plot.GetRight() - width - gap, plot.Y + gap, width, height},                     // top right
            {plot.GetRight() - width - gap, plot.GetBottom() - height - gap, width, height},  // bottom right
            {plot.X + gap, plot.Y + gap, width, height},                                      // top left
            {plot.X + gap, plot.GetBottom() - height - gap, width, height}}};                 // bottom left
    };
    std::array<Gdiplus::RectF, 4> legendCorners{}, insetCorners{};
    if (legendSize) legendCorners = corners(legendSize->Width, legendSize->Height, 6.0f);
    if (insetSize) insetCorners = corners(insetSize->Width, insetSize->Height, 8.0f);
    // Important points, in pixels.
    std::vector<Gdiplus::PointF> important;
    const double xSpan = view.xMax - view.xMin;
    const double yLow = std::min(view.yMin, view.yMax), yHigh = std::max(view.yMin, view.yMax);
    auto pixel = [&](double x, double y) {
        return Gdiplus::PointF(static_cast<float>(view.left + view.width * (x - view.xMin) / xSpan),
                               static_cast<float>(plotPixelY(view, y)));
    };
    // The figure's highest and lowest points matter only when they are prominent: clear of
    // the extremes reached elsewhere (repetitive peaks and troughs are not a reason to move).
    constexpr int bins = 24;
    struct Extreme { double y; double x; bool set = false; };
    std::array<Extreme, bins> binMax{}, binMin{};
    for (const Trace& trace : figure.traces) {
        if (trace.heatmap || !(xSpan > 0)) continue;
        const std::size_t count = std::min(trace.x.size(), trace.y.size());
        for (std::size_t index = 0; index < count; ++index) {
            const double x = trace.x[index], y = trace.y[index];
            if (!std::isfinite(x) || !std::isfinite(y) || x < view.xMin || x > view.xMax || y < yLow || y > yHigh) continue;
            const int bin = std::clamp(static_cast<int>((x - view.xMin) / xSpan * bins), 0, bins - 1);
            if (!binMax[bin].set || y > binMax[bin].y) binMax[bin] = {y, x, true};
            if (!binMin[bin].set || y < binMin[bin].y) binMin[bin] = {y, x, true};
            if (detail && detail->transientWindow && x >= detail->xMin && x <= detail->xMax)
                important.push_back(pixel(x, y));
        }
    }
    auto prominent = [&](const std::array<Extreme, bins>& extremes, bool highest) {
        int best = -1;
        for (int bin = 0; bin < bins; ++bin)
            if (extremes[bin].set && (best < 0 || (highest ? extremes[bin].y > extremes[best].y : extremes[bin].y < extremes[best].y)))
                best = bin;
        if (best < 0) return;
        double runnerUp = std::numeric_limits<double>::quiet_NaN();
        for (int bin = 0; bin < bins; ++bin) {
            if (!extremes[bin].set || std::abs(bin - best) <= 1) continue;
            if (std::isnan(runnerUp) || (highest ? extremes[bin].y > runnerUp : extremes[bin].y < runnerUp)) runnerUp = extremes[bin].y;
        }
        if (std::isnan(runnerUp) || std::abs(extremes[best].y - runnerUp) > 0.05 * (yHigh - yLow))
            important.push_back(pixel(extremes[best].x, extremes[best].y));
    };
    prominent(binMax, true); prominent(binMin, false);
    auto hidden = [&](const Gdiplus::RectF& area) {
        double n = 0;
        for (const auto& point : important) if (area.Contains(point)) n += 1;
        return n;
    };
    OverlayLayout result;
    double best = std::numeric_limits<double>::infinity();
    for (int l = 0; l < (legendSize ? 4 : 1); ++l)
        for (int i = 0; i < (insetSize ? 4 : 1); ++i) {
            if (legendSize && insetSize && l == i) continue;
            double hides = 0, moves = 0;
            if (legendSize) { hides += hidden(legendCorners[l]); moves += l != 0; }
            if (insetSize) { hides += hidden(insetCorners[i]); moves += i != 1; }
            const double cost = hides * 1000.0 + moves;
            if (cost < best) {
                best = cost; result.cost = cost;
                result.legend.reset(); result.inset.reset();
                if (legendSize) result.legend = legendCorners[l];
                if (insetSize) result.inset = insetCorners[i];
            }
        }
    return result;
}

std::optional<DetailInset> resolvedInset(const Figure& figure) {
    return figure.inset ? figure.inset : ordinateScale(figure).kind == ValueScaleKind::linear
        ? automaticDetailInset(figure) : std::nullopt;
}

void inset(Gdiplus::Graphics& graphics, const Figure& figure, const std::optional<DetailInset>& resolved,
           const Gdiplus::RectF& area, const theme::Fonts& fonts) {
    if (!resolved) return;
    // Opaque: main-panel curves must not show through and read as inset data.
    Gdiplus::SolidBrush background(theme::color(theme::plotBackground));
    Gdiplus::Pen outline(theme::color(theme::accent, 190), 1.0f);
    graphics.FillRectangle(&background, area);
    graphics.DrawRectangle(&outline, area);
    // The inset is a small plot in its own right: title, gridded axes with tick labels
    // and axis titles, so its window and scale are unambiguous.
    const float s = scientificTextScale(fonts.body);
    const float left = 80 * s, bottom = 36 * s, top = 20 * s, right = 10 * s;
    const Gdiplus::RectF inner(area.X + left, area.Y + top,
        std::max(1.0f, area.Width - left - right), std::max(1.0f, area.Height - top - bottom));
    Range xr{resolved->xMin, resolved->xMax};
    Range yr{resolved->yMin, resolved->yMax};
    std::wstring heading = resolved->transientWindow ? L"detail: start-up transient"
        : resolved->automatic ? L"automatic detail" : L"detail";
    if (resolved->symmetricLog) {
        wchar_t core[64]{};
        swprintf_s(core, L"; symmetric log scale, linear within \u00b1%.3g", resolved->linearThreshold);
        heading += core;
    }
    drawText(graphics, heading, {area.X + 7, area.Y, area.Width - 14, top},
             fonts.caption, theme::accentStrong);
    auto exponent = [](const Range& range) {
        const auto magnitude = std::max(std::abs(range.low), std::abs(range.high));
        return magnitude > 0 && (magnitude >= 10000 || magnitude < .001)
            ? static_cast<int>(std::floor(std::log10(magnitude))) : 0;
    };
    ValueScale yScale = ordinateScale(figure);
    if (resolved->symmetricLog)
        yScale = ValueScale{ValueScaleKind::symmetric_log, resolved->linearThreshold, "inset", resolved->decadeWeight};
    const bool linearY = yScale.kind == ValueScaleKind::linear;
    const int xExponent = exponent(xr), yExponent = linearY ? exponent(yr) : 0;
    Gdiplus::Pen gridPen(theme::color(theme::plotGrid), 1.0f);
    Gdiplus::Pen axisPen(theme::color(theme::plotMuted, 200), 1.0f);
    wchar_t label[48]{};
    const PlotTransform insetView{inner.X, inner.Y, inner.Width, inner.Height, xr.low, xr.high, yr.low, yr.high, yScale};
    constexpr int divisions = 3;
    for (int index = 0; index <= divisions; ++index) {
        const float fraction = static_cast<float>(index) / divisions;
        const float gx = inner.X + fraction * inner.Width;
        graphics.DrawLine(&gridPen, gx, inner.Y, gx, inner.GetBottom());
        graphics.DrawLine(&axisPen, gx, inner.GetBottom(), gx, inner.GetBottom() + 3 * s);
        swprintf_s(label, L"%.3g", (xr.low + fraction * (xr.high - xr.low)) / std::pow(10., xExponent));
        drawMathText(graphics, label, {gx - 30 * s, inner.GetBottom() + 3 * s, 60 * s, 15 * s},
            fonts.caption, theme::plotMuted, Gdiplus::StringAlignmentCenter);
    }
    // Ordinate ticks: evenly spaced on a linear scale; 0 and signed decades on symmetric log.
    std::vector<std::pair<double, std::wstring>> yTicks;
    if (linearY) {
        for (int index = 0; index <= divisions; ++index) {
            const double value = yr.high - (yr.high - yr.low) * index / divisions;
            swprintf_s(label, L"%.3g", value / std::pow(10., yExponent));
            yTicks.emplace_back(value, label);
        }
    } else {
        // 0, the edges of the linear core (plain values), then decades clear of the core.
        const double reach = std::max(std::abs(yr.low), std::abs(yr.high));
        const double core = yScale.linearThreshold;
        if (yr.low <= 0 && yr.high >= 0) yTicks.emplace_back(0.0, L"0");
        swprintf_s(label, L"%.3g", core);
        if (core <= yr.high) yTicks.emplace_back(core, label);
        swprintf_s(label, L"-%.3g", core);
        if (-core >= yr.low) yTicks.emplace_back(-core, label);
        const int first = static_cast<int>(std::ceil(std::log10(std::max(5.0 * core, 1e-300))));
        const int last = static_cast<int>(std::floor(std::log10(std::max(reach, 1e-300))));
        std::vector<int> decades;
        for (int k = first; k <= last; ++k) decades.push_back(k);
        const std::size_t stride = decades.size() > 2 ? (decades.size() + 1) / 2 : 1;
        for (std::size_t n = decades.size(); n-- > 0;) {
            if ((decades.size() - 1 - n) % stride) continue;
            const int k = decades[n];
            const double magnitude = std::pow(10., k);
            const std::wstring power = L"10^{" + std::to_wstring(k) + L"}";
            if (magnitude <= yr.high) yTicks.emplace_back(magnitude, power);
            if (-magnitude >= yr.low) yTicks.emplace_back(-magnitude, L"-" + power);
        }
    }
    std::vector<float> placed;  // ticks in priority order; skip any that would crowd a placed one
    for (const auto& [value, text] : yTicks) {
        const float gy = static_cast<float>(plotPixelY(insetView, value));
        if (gy < inner.Y - 1 || gy > inner.GetBottom() + 1) continue;
        if (std::any_of(placed.begin(), placed.end(), [&](float other) { return std::abs(other - gy) < 14 * s; })) continue;
        placed.push_back(gy);
        graphics.DrawLine(&gridPen, inner.X, gy, inner.GetRight(), gy);
        graphics.DrawLine(&axisPen, inner.X - 3 * s, gy, inner.X, gy);
        drawMathText(graphics, text, {inner.X - 62 * s, gy - 8 * s, 58 * s, 16 * s},
            fonts.caption, theme::plotMuted, Gdiplus::StringAlignmentFar);
    }
    graphics.DrawRectangle(&axisPen, inner);
    auto axisTitle = [](const std::string& name, int power) {
        return wide(scientificMathLabel(name) + (power
            ? "\\; \\times 10^{" + std::to_string(power) + "}" : ""));
    };
    drawMathText(graphics, axisTitle(figure.xLabel, xExponent),
        {inner.X, inner.GetBottom() + 17 * s, inner.Width, 17 * s},
        fonts.caption, theme::plotMuted, Gdiplus::StringAlignmentCenter);
    drawVerticalText(graphics, axisTitle(figure.yLabel, yExponent),
        {area.X + 2 * s, inner.Y, 15 * s, inner.Height}, fonts.caption, theme::plotMuted);
    const Gdiplus::GraphicsState state = graphics.Save();
    graphics.SetClip(inner, Gdiplus::CombineModeIntersect);
    for (std::size_t index = 0; index < figure.traces.size(); ++index) {
        if (!resolved->traceIndices.empty() && std::find(resolved->traceIndices.begin(),
                resolved->traceIndices.end(), index) == resolved->traceIndices.end()) continue;
        drawPlotSeries(graphics, figure.traces[index], index,
            insetView);
    }
    graphics.Restore(state);
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
        // Non-trajectory 3-D plots (stem plots) with several labelled series get a
        // legend. The camera controls occupy the top right and the interaction hint the
        // bottom edge, so it sits bottom right, just above the hint.
        if (figure.boxAxes3D) {
            std::size_t labelled = 0;
            for (const auto& trace : figure.traces) if (!trace.label.empty()) ++labelled;
            if (labelled >= 2) {
                const auto size = plotLegendSize(graphics, figure, plot, fonts);
                const Gdiplus::PointF origin(plot.X + plot.Width - size.Width - 8.0f,
                                             plot.Y + plot.Height - size.Height - 34.0f);
                drawPlotLegend(graphics, figure, plot, fonts, &origin);
            }
        }
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
    const bool ordinaryLegend = figure.kind != "krom_spectrum" && figure.kind != "stem"
        && !(figure.kind == "krom_time" && figure.traces.size() == 1)
        && !figure.traces.empty() && !figure.traces.front().heatmap;
    Figure selected = figure;
    if (figure.kind == "krom_spectrum") {
        selected.traces.clear();
        for (const auto& t : figure.traces) if (t.selected) selected.traces.push_back(t);
    }
    const Figure* legendFigure = ordinaryLegend ? &figure : figure.kind == "krom_spectrum" ? &selected : nullptr;
    std::optional<Gdiplus::SizeF> legendSize;
    if (legendFigure) {
        const auto size = plotLegendSize(graphics, *legendFigure, plot, fonts);
        if (size.Width > 0 && size.Height > 0) legendSize = size;
    }
    const auto detail = resolvedInset(figure);
    // Full-size inset unless only a smaller one avoids hiding something that matters.
    OverlayLayout layout = overlayLayout(figure, plot, transform, legendSize, std::nullopt, detail);
    if (detail) {
        double best = std::numeric_limits<double>::infinity();
        for (const float factor : {1.0f, 0.82f, 0.68f}) {
            const auto candidate = overlayLayout(figure, plot, transform, legendSize,
                Gdiplus::SizeF(plot.Width * 0.40f * factor, plot.Height * 0.42f * factor), detail);
            const double cost = candidate.cost + 10.0 * (1.0 - factor);
            if (cost < best) { best = cost; layout = candidate; }
        }
    }
    if (legendFigure && layout.legend) {
        const Gdiplus::PointF origin(layout.legend->X, layout.legend->Y);
        drawPlotLegend(graphics, *legendFigure, plot, fonts, &origin);
    }
    if (layout.inset) inset(graphics, figure, detail, *layout.inset, fonts);
}

}
