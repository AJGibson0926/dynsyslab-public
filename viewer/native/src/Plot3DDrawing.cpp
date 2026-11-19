#include "Plot3DDrawing.hpp"

#include "MathTextDrawing.hpp"
#include "PlotFigureDrawing.hpp"
#include "PlotSeries.hpp"
#include "dynsyslab_viewer/MathText.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <vector>

namespace dsl {
namespace {

struct Segment {
    Plot3DProjectedPoint first;
    Plot3DProjectedPoint second;
    std::size_t traceIndex{};
    std::size_t sampleIndex{};
    double depth{};
};

Gdiplus::PointF point(Plot3DProjectedPoint value) {
    return {static_cast<float>(value.x), static_cast<float>(value.y)};
}

std::pair<double, double> finiteRange(const std::vector<double>& values) {
    double low = std::numeric_limits<double>::infinity();
    double high = -low;
    for (const auto value : values) if (std::isfinite(value)) {
        low = std::min(low, value); high = std::max(high, value);
    }
    return std::isfinite(low) ? std::make_pair(low, high)
                              : std::make_pair(0.0, 1.0);
}

void drawArrow(Gdiplus::Graphics& graphics, const Plot3DControlBounds& area,
               Plot3DControl control, COLORREF color) {
    const float cx = static_cast<float>(area.left + area.width * 0.5);
    const float cy = static_cast<float>(area.top + area.height * 0.5);
    float dx = 0.0f, dy = 0.0f;
    if (control == Plot3DControl::orbit_left || control == Plot3DControl::pan_left) dx = -1.0f;
    if (control == Plot3DControl::orbit_right || control == Plot3DControl::pan_right) dx = 1.0f;
    if (control == Plot3DControl::orbit_up || control == Plot3DControl::pan_up) dy = -1.0f;
    if (control == Plot3DControl::orbit_down || control == Plot3DControl::pan_down) dy = 1.0f;
    Gdiplus::Pen pen(theme::color(color), 1.8f);
    pen.SetStartCap(Gdiplus::LineCapRound); pen.SetEndCap(Gdiplus::LineCapRound);
    const Gdiplus::PointF start(cx - dx * 8.0f, cy - dy * 8.0f);
    const Gdiplus::PointF end(cx + dx * 8.0f, cy + dy * 8.0f);
    graphics.DrawLine(&pen, start, end);
    const float px = -dy, py = dx;
    graphics.DrawLine(&pen, end, {end.X - dx * 6.0f + px * 4.5f,
                                  end.Y - dy * 6.0f + py * 4.5f});
    graphics.DrawLine(&pen, end, {end.X - dx * 6.0f - px * 4.5f,
                                  end.Y - dy * 6.0f - py * 4.5f});
    const bool orbit = control == Plot3DControl::orbit_left
        || control == Plot3DControl::orbit_right
        || control == Plot3DControl::orbit_up
        || control == Plot3DControl::orbit_down;
    if (orbit) graphics.DrawArc(&pen, cx - 12.0f, cy - 12.0f, 24.0f, 24.0f,
                                dx < 0 || dy < 0 ? 35.0f : 215.0f, 110.0f);
}

void drawProjectionIcon(Gdiplus::Graphics& graphics,
                        const Plot3DControlBounds& area, COLORREF color) {
    Gdiplus::Pen pen(theme::color(color), 1.5f);
    const float x = static_cast<float>(area.left + 11.0);
    const float y = static_cast<float>(area.top + 11.0);
    graphics.DrawRectangle(&pen, x, y + 5.0f, 14.0f, 14.0f);
    graphics.DrawRectangle(&pen, x + 6.0f, y, 14.0f, 14.0f);
    graphics.DrawLine(&pen, x, y + 5.0f, x + 6.0f, y);
    graphics.DrawLine(&pen, x + 14.0f, y + 5.0f, x + 20.0f, y);
    graphics.DrawLine(&pen, x + 14.0f, y + 19.0f, x + 20.0f, y + 14.0f);
}

void drawHomeIcon(Gdiplus::Graphics& graphics,
                  const Plot3DControlBounds& area, COLORREF color) {
    Gdiplus::Pen pen(theme::color(color), 1.7f);
    const float cx = static_cast<float>(area.left + area.width * 0.5);
    const float cy = static_cast<float>(area.top + area.height * 0.5);
    Gdiplus::PointF roof[3]{{cx - 10, cy - 1}, {cx, cy - 10}, {cx + 10, cy - 1}};
    graphics.DrawLines(&pen, roof, 3);
    graphics.DrawRectangle(&pen, cx - 7.0f, cy - 1.0f, 14.0f, 11.0f);
}

void drawGizmos(Gdiplus::Graphics& graphics, const Plot3DViewport& viewport,
                const Plot3DCamera& camera, Plot3DControl hovered,
                const theme::Fonts& fonts) {
    const auto controls = plot3DControls(viewport);
    for (const auto& control : controls) {
        const bool active = control.control == hovered;
        Gdiplus::SolidBrush fill(theme::color(active ? theme::selection
                                                     : theme::panelRaised, active ? 250 : 224));
        Gdiplus::Pen border(theme::color(active ? theme::focus : theme::plotBorder,
                                         active ? 255 : 210), active ? 1.8f : 1.0f);
        const Gdiplus::RectF area(static_cast<float>(control.left),
            static_cast<float>(control.top), static_cast<float>(control.width),
            static_cast<float>(control.height));
        graphics.FillRectangle(&fill, area); graphics.DrawRectangle(&border, area);
        const COLORREF icon = active ? theme::accentStrong : theme::plotMuted;
        if (control.control >= Plot3DControl::orbit_left
            && control.control <= Plot3DControl::pan_right)
            drawArrow(graphics, control, control.control, icon);
        else if (control.control == Plot3DControl::zoom_in
                 || control.control == Plot3DControl::zoom_out) {
            Gdiplus::Pen pen(theme::color(icon), 1.8f);
            const float cx = area.X + area.Width * 0.5f;
            const float cy = area.Y + area.Height * 0.5f;
            graphics.DrawLine(&pen, cx - 8, cy, cx + 8, cy);
            if (control.control == Plot3DControl::zoom_in)
                graphics.DrawLine(&pen, cx, cy - 8, cx, cy + 8);
        } else if (control.control == Plot3DControl::toggle_projection)
            drawProjectionIcon(graphics, control, icon);
        else if (control.control == Plot3DControl::reset)
            drawHomeIcon(graphics, control, icon);
    }
    if (hovered != Plot3DControl::none) {
        const auto found = std::find_if(controls.begin(), controls.end(),
            [&](const auto& value) { return value.control == hovered; });
        if (found != controls.end()) {
            const auto label = wide(plot3DControlName(hovered));
            const float width = 132.0f, height = 25.0f;
            const float left = std::clamp(static_cast<float>(found->left + found->width - width),
                static_cast<float>(viewport.left),
                static_cast<float>(viewport.left + viewport.width - width));
            const float top = static_cast<float>(found->top + found->height + 4.0);
            Gdiplus::SolidBrush fill(theme::color(theme::text, 242));
            graphics.FillRectangle(&fill, left, top, width, height);
            drawText(graphics, label, {left + 5, top, width - 10, height},
                     fonts.caption, RGB(255, 255, 255), Gdiplus::StringAlignmentCenter);
        }
    }
    (void)camera;
}

void drawTriad(Gdiplus::Graphics& graphics, const Plot3DCamera& camera,
               const Plot3DViewport& viewport, const theme::Fonts& fonts) {
    auto triadState = camera.state();
    triadState.panX = 0.0; triadState.panY = 0.0; triadState.zoom = 0.72;
    const Plot3DCamera triadCamera(triadState);
    const Plot3DBounds unit{{-1, -1, -1}, {1, 1, 1}};
    const Plot3DViewport area{viewport.left + 12.0,
        viewport.top + viewport.height - 88.0, 76.0, 76.0};
    const auto origin = triadCamera.project({0, 0, 0}, unit, area);
    const std::array<Plot3DPoint, 3> axes{{{0.75, 0, 0}, {0, 0.75, 0}, {0, 0, 0.75}}};
    const std::array<COLORREF, 3> colors{{theme::accentStrong, theme::blue, theme::danger}};
    const wchar_t* labels[3] = {L"x", L"y", L"z"};
    for (std::size_t index = 0; index < axes.size(); ++index) {
        const auto end = triadCamera.project(axes[index], unit, area);
        Gdiplus::Pen pen(theme::color(colors[index]), 2.2f);
        graphics.DrawLine(&pen, point(origin), point(end));
        drawMathText(graphics, labels[index],
            {static_cast<float>(end.x - 8), static_cast<float>(end.y - 10), 18, 18},
            fonts.caption, colors[index], Gdiplus::StringAlignmentCenter);
    }
}

std::wstring cameraText(const Plot3DCameraState& state) {
    std::wostringstream value;
    value << std::fixed << std::setprecision(0) << L"Az " << state.azimuthDegrees
          << L"\u00b0  El " << state.elevationDegrees << L"\u00b0  ";
    value << std::setprecision(2) << L"Zoom " << state.zoom << L"\u00d7  ";
    value << (state.projection == Plot3DProjection::perspective
        ? L"Perspective" : L"Orthographic");
    return value.str();
}

}  // namespace

// About `target` round-valued ticks (1, 2 or 5 times a power of ten) inside [low, high].
static std::vector<double> niceTicks(double low, double high, int target) {
    std::vector<double> ticks;
    if (!(high > low) || !std::isfinite(low) || !std::isfinite(high)) return ticks;
    const double raw = (high - low) / std::max(1, target);
    const double power = std::pow(10.0, std::floor(std::log10(raw)));
    double step = power;
    for (const double m : {1.0, 2.0, 5.0, 10.0}) { step = m * power; if (step >= raw) break; }
    for (double v = std::ceil(low / step) * step; v <= high + 1e-9 * step; v += step)
        ticks.push_back(std::abs(v) < 1e-12 * step ? 0.0 : v);
    return ticks;
}

static std::wstring tickText(double value) {
    wchar_t text[40]{};
    const double magnitude = std::abs(value);
    if (magnitude != 0 && (magnitude >= 1e4 || magnitude < 1e-3)) {
        const int k = static_cast<int>(std::floor(std::log10(magnitude)));
        const double mantissa = value / std::pow(10.0, k);
        if (std::abs(std::abs(mantissa) - 1.0) < 1e-9) return (value < 0 ? L"-10^{" : L"10^{") + std::to_wstring(k) + L"}";
        swprintf_s(text, L"%.2g\times10^{%d}", mantissa, k);
        return text;
    }
    swprintf_s(text, L"%.4g", value);
    return text;
}

// Centred text on a near-opaque plot-coloured backing, so labels stay readable where
// they cross stems or curves.
static void drawBackedText(Gdiplus::Graphics& graphics, const std::wstring& text, const Gdiplus::RectF& box,
                           HFONT font, COLORREF color) {
    const auto size = measureMathText(graphics, text, font);
    const COLORREF fill = theme::plotBackground;
    Gdiplus::SolidBrush backing(Gdiplus::Color(225, GetRValue(fill), GetGValue(fill), GetBValue(fill)));
    const float w = std::min(box.Width, size.Width + 6.0f), h = std::min(box.Height, size.Height + 2.0f);
    graphics.FillRectangle(&backing, box.X + (box.Width - w) / 2, box.Y + (box.Height - h) / 2, w, h);
    drawMathText(graphics, text, box, font, color, Gdiplus::StringAlignmentCenter);
}

// Tick marks and labels along one box axis through `start` (data coordinates; the
// `axis` component varies). Labels sit on the side away from the box centre; height
// ticks point sideways. Ticks within 22 px of `skipNear` (a corner another axis
// already labels) are left out.
static void drawAxisTicks(Gdiplus::Graphics& graphics, const Plot3DCamera& camera, const Plot3DBounds& bounds,
                   const Plot3DViewport& viewport, Plot3DPoint start, int axis,
                   const std::vector<std::pair<double, std::wstring>>& ticks, COLORREF color,
                   const theme::Fonts& fonts, const Plot3DPoint* skipNear = nullptr) {
    const Plot3DPoint centre{(bounds.minimum.x + bounds.maximum.x) / 2, (bounds.minimum.y + bounds.maximum.y) / 2,
                             (bounds.minimum.z + bounds.maximum.z) / 2};
    const auto mid = camera.project(centre, bounds, viewport);
    Gdiplus::Pen pen(theme::color(color), 1.2f);
    const auto origin = camera.project(skipNear ? *skipNear : start, bounds, viewport);
    float lastX = -1e9f, lastY = -1e9f;
    for (const auto& [value, text] : ticks) {
        Plot3DPoint p = start;
        (axis == 0 ? p.x : axis == 1 ? p.y : p.z) = value;
        const auto at = camera.project(p, bounds, viewport);
        if (!at.visible) continue;
        if (skipNear && std::hypot(at.x - origin.x, at.y - origin.y) < 22) continue;
        double dx = at.x - mid.x, dy = at.y - mid.y;
        if (axis == 2) { dx = at.x <= mid.x ? -1 : 1; dy = 0; }   // height ticks point sideways
        const double length = std::max(1e-9, std::hypot(dx, dy));
        dx /= length; dy /= length;
        const float tx = static_cast<float>(at.x + 6 * dx), ty = static_cast<float>(at.y + 6 * dy);
        if (std::hypot(tx - lastX, ty - lastY) < 26) continue;   // keep labels apart
        lastX = tx; lastY = ty;
        graphics.DrawLine(&pen, static_cast<float>(at.x), static_cast<float>(at.y), tx, ty);
        const float lx = static_cast<float>(at.x + 26 * dx), ly = static_cast<float>(at.y + 16 * dy);
        drawBackedText(graphics, text, {lx - 40, ly - 10, 80, 20}, fonts.caption, theme::plotMuted);
    }
}

void drawPlot3D(Gdiplus::Graphics& graphics, const Figure& sourceFigure,
                const Plot3DViewport& viewport, const Plot3DCamera& camera,
                Plot3DControl hoveredControl, const theme::Fonts& fonts) {
    // Box-axis plots (stems) honour the value-scale control on their height axis.
    const auto heightScale = plot3DHeightScale(sourceFigure);
    const Figure figure = plot3DDisplayFigure(sourceFigure);
    const auto bounds = plot3DBounds(figure);
    const auto saved = graphics.Save();
    const Gdiplus::RectF clip(static_cast<float>(viewport.left),
        static_cast<float>(viewport.top), static_cast<float>(viewport.width),
        static_cast<float>(viewport.height));
    graphics.SetClip(clip, Gdiplus::CombineModeIntersect);
    Gdiplus::SolidBrush plotFill(theme::color(theme::plotBackground));
    graphics.FillRectangle(&plotFill, clip);

    const std::array<Plot3DPoint, 8> corners{{
        {bounds.minimum.x, bounds.minimum.y, bounds.minimum.z},
        {bounds.maximum.x, bounds.minimum.y, bounds.minimum.z},
        {bounds.minimum.x, bounds.maximum.y, bounds.minimum.z},
        {bounds.maximum.x, bounds.maximum.y, bounds.minimum.z},
        {bounds.minimum.x, bounds.minimum.y, bounds.maximum.z},
        {bounds.maximum.x, bounds.minimum.y, bounds.maximum.z},
        {bounds.minimum.x, bounds.maximum.y, bounds.maximum.z},
        {bounds.maximum.x, bounds.maximum.y, bounds.maximum.z}}};
    std::array<Plot3DProjectedPoint, 8> projected{};
    for (std::size_t index = 0; index < corners.size(); ++index)
        projected[index] = camera.project(corners[index], bounds, viewport);
    constexpr int edges[12][2] = {{0,1},{0,2},{1,3},{2,3},{4,5},{4,6},
        {5,7},{6,7},{0,4},{1,5},{2,6},{3,7}};
    Gdiplus::Pen boxPen(theme::color(theme::plotGrid), 1.0f);
    for (const auto& edge : edges)
        graphics.DrawLine(&boxPen, point(projected[edge[0]]), point(projected[edge[1]]));

    std::vector<Segment> segments;
    for (std::size_t traceIndex = 0; traceIndex < figure.traces.size(); ++traceIndex) {
        const auto& trace = figure.traces[traceIndex];
        const auto count = std::min({trace.x.size(), trace.y.size(), trace.z.size()});
        for (std::size_t index = 1; index < count; ++index) {
            const Plot3DPoint firstValue{trace.x[index - 1], trace.y[index - 1], trace.z[index - 1]};
            const Plot3DPoint secondValue{trace.x[index], trace.y[index], trace.z[index]};
            const auto first = camera.project(firstValue, bounds, viewport);
            const auto second = camera.project(secondValue, bounds, viewport);
            if (first.visible && second.visible)
                segments.push_back({first, second, traceIndex, index,
                                    (first.depth + second.depth) * 0.5});
        }
    }
    std::stable_sort(segments.begin(), segments.end(),
        [](const auto& first, const auto& second) { return first.depth < second.depth; });
    for (const auto& segment : segments) {
        const auto& trace = figure.traces[segment.traceIndex];
        Gdiplus::Color color = plotSeriesColor(trace, segment.traceIndex);
        if (trace.encoding == "time" && !trace.colorValues.empty()) {
            const auto range = finiteRange(trace.colorValues);
            const auto sample = std::min(segment.sampleIndex - 1, trace.colorValues.size() - 1);
            const double value = trace.colorValues[sample];
            color = applyPlotSeriesOpacity(trajectoryTimeColor(range.second > range.first
                ? (value - range.first) / (range.second - range.first) : 0.0), trace);
        }
        const double normalizedDepth = std::clamp((segment.depth + 1.75) / 3.5, 0.0, 1.0);
        color = Gdiplus::Color(static_cast<BYTE>(color.GetA() *
            (150.0 + normalizedDepth * 105.0) / 255.0),
            color.GetR(), color.GetG(), color.GetB());
        Gdiplus::Pen pen(color);
        stylePlotSeriesPen(pen, trace);
        graphics.DrawLine(&pen, point(segment.first), point(segment.second));
    }

    const std::array<COLORREF, 3> axisColors{{theme::accentStrong, theme::blue, theme::danger}};
    const std::array<std::string, 3> axisLabels{{figure.xLabel, figure.yLabel,
        (figure.zLabel.empty() ? std::string("z") : figure.zLabel)
            + (heightScale.logarithmic ? std::string(" (log)") : std::string())}};
    // Axis lines (start, end). Trajectory plots keep one shared corner. Box plots (stems)
    // put each base axis on the floor edge nearest the viewer, so stems never hide it,
    // and the height axis on the floor corner furthest left on screen.
    std::array<Plot3DPoint, 3> axisStarts{}, axisEnds{};
    const double x0 = bounds.minimum.x, x1 = bounds.maximum.x, y0 = bounds.minimum.y, y1 = bounds.maximum.y;
    const double z0 = bounds.minimum.z, z1 = bounds.maximum.z;
    Plot3DPoint baseCorner{x0, y0, z0};
    if (figure.boxAxes3D) {
        const auto depth = [&](Plot3DPoint p) { return camera.project(p, bounds, viewport).depth; };
        const double yFront = depth({(x0 + x1) / 2, y0, z0}) >= depth({(x0 + x1) / 2, y1, z0}) ? y0 : y1;
        const double xFront = depth({x0, (y0 + y1) / 2, z0}) >= depth({x1, (y0 + y1) / 2, z0}) ? x0 : x1;
        axisStarts[0] = {x0, yFront, z0}; axisEnds[0] = {x1, yFront, z0};
        axisStarts[1] = {xFront, y0, z0}; axisEnds[1] = {xFront, y1, z0};
        baseCorner = {xFront, yFront, z0};
        Plot3DPoint left{x0, y0, z0}; double leftX = std::numeric_limits<double>::infinity();
        for (const Plot3DPoint c : {Plot3DPoint{x0, y0, z0}, Plot3DPoint{x1, y0, z0}, Plot3DPoint{x0, y1, z0}, Plot3DPoint{x1, y1, z0}}) {
            const double sx = camera.project(c, bounds, viewport).x;
            if (sx < leftX) { leftX = sx; left = c; }
        }
        axisStarts[2] = left; axisEnds[2] = {left.x, left.y, z1};
    } else {
        axisStarts = {Plot3DPoint{x0, y0, z0}, Plot3DPoint{x0, y0, z0}, Plot3DPoint{x0, y0, z0}};
        axisEnds = {Plot3DPoint{x1, y0, z0}, Plot3DPoint{x0, y1, z0}, Plot3DPoint{x0, y0, z1}};
    }
    for (std::size_t index = 0; index < 3; ++index) {
        const auto start = camera.project(axisStarts[index], bounds, viewport);
        const auto end = camera.project(axisEnds[index], bounds, viewport);
        Gdiplus::Pen pen(theme::color(axisColors[index]), 2.0f);
        graphics.DrawLine(&pen, point(start), point(end));
        const Gdiplus::RectF titleBox{static_cast<float>(end.x - 60), static_cast<float>(end.y - 25), 120, 22};
        if (figure.boxAxes3D) drawBackedText(graphics, wide(scientificMathLabel(axisLabels[index])), titleBox, fonts.body, axisColors[index]);
        else drawMathText(graphics, wide(scientificMathLabel(axisLabels[index])), titleBox,
            fonts.body, axisColors[index], Gdiplus::StringAlignmentCenter);
    }
    if (figure.boxAxes3D) {   // tick marks with physical values on all three axes
        const double lows[3] = {bounds.minimum.x, bounds.minimum.y, bounds.minimum.z};
        const double highs[3] = {bounds.maximum.x, bounds.maximum.y, bounds.maximum.z};
        for (int axis = 0; axis < 3; ++axis) {
            std::vector<std::pair<double, std::wstring>> ticks;
            if (axis == 2 && heightScale.logarithmic) {
                const int first = static_cast<int>(std::ceil(lows[2] - 1e-9));
                const int last = static_cast<int>(std::floor(highs[2] + 1e-9));
                const int stride = std::max(1, (last - first + 5) / 6);
                for (int k = last; k >= first; k -= stride)
                    ticks.emplace_back(double(k), L"10^{" + std::to_wstring(k) + L"}");
            } else {
                for (const double v : niceTicks(lows[axis], highs[axis], 4)) ticks.emplace_back(v, tickText(v));
            }
            // The x axis labels the corner it shares with the y axis; the height axis
            // leaves out its floor tick (a base axis may label that corner).
            const Plot3DPoint* skip = axis == 1 ? &baseCorner : axis == 2 ? &axisStarts[2] : nullptr;
            drawAxisTicks(graphics, camera, bounds, viewport, axisStarts[axis], axis, ticks, axisColors[axis], fonts, skip);
        }
    }

    (void)sourceFigure;
    if (!figure.traces.empty() && !figure.boxAxes3D) {  // IC marker: trajectories only
        const auto& trace = figure.traces.front();
        if (!trace.x.empty() && !trace.y.empty() && !trace.z.empty()) {
            const auto initial = camera.project({trace.x.front(), trace.y.front(), trace.z.front()},
                                                bounds, viewport);
            Gdiplus::SolidBrush white(theme::color(theme::plotBackground));
            Gdiplus::Pen marker(theme::color(theme::danger), 1.6f);
            graphics.FillEllipse(&white, static_cast<float>(initial.x - 5),
                                 static_cast<float>(initial.y - 5), 10.0f, 10.0f);
            graphics.DrawEllipse(&marker, static_cast<float>(initial.x - 5),
                                 static_cast<float>(initial.y - 5), 10.0f, 10.0f);
            drawText(graphics, L"IC", {static_cast<float>(initial.x + 7),
                static_cast<float>(initial.y - 10), 28, 20}, fonts.caption,
                theme::plotMuted);
        }
    }

    drawTriad(graphics, camera, viewport, fonts);
    drawText(graphics, cameraText(camera.state()),
        {static_cast<float>(viewport.left + 10), static_cast<float>(viewport.top + 7),
         static_cast<float>(std::max(1.0, viewport.width - 204)), 24},
        fonts.caption, theme::plotMuted);
    drawText(graphics,
        L"Drag: rotate   Ctrl / middle-drag: move   Alt / right-drag: zoom   Wheel: zoom   Home: reset",
        {static_cast<float>(viewport.left + 10),
         static_cast<float>(viewport.top + viewport.height - 25),
         static_cast<float>(std::max(1.0, viewport.width - 20)), 20},
        fonts.caption, theme::plotMuted, Gdiplus::StringAlignmentCenter);
    drawGizmos(graphics, viewport, camera, hoveredControl, fonts);
    graphics.Restore(saved);
}

}  // namespace dsl
