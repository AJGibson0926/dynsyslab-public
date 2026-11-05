#include "dynsyslab_viewer/Plot3D.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace dsl {
namespace {

constexpr double pi = 3.14159265358979323846;

double radians(double degrees) { return degrees * pi / 180.0; }

double finiteOr(double value, double fallback) {
    return std::isfinite(value) ? value : fallback;
}

double dot(Plot3DPoint a, Plot3DPoint b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Plot3DPoint subtract(Plot3DPoint a, Plot3DPoint b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

Plot3DPoint interpolate(Plot3DPoint a, Plot3DPoint b, double fraction) {
    return {a.x + fraction * (b.x - a.x),
            a.y + fraction * (b.y - a.y),
            a.z + fraction * (b.z - a.z)};
}

Plot3DPoint tracePoint(const Trace& trace, std::size_t index) {
    return {trace.x[index], trace.y[index], trace.z[index]};
}

bool finite(Plot3DPoint point) {
    return std::isfinite(point.x) && std::isfinite(point.y)
        && std::isfinite(point.z);
}

Plot3DProbe segmentProbe(const Plot3DPoint& firstValue,
                         const Plot3DPoint& secondValue,
                         Plot3DProjectedPoint first,
                         Plot3DProjectedPoint second,
                         std::size_t traceIndex,
                         std::size_t segmentIndex,
                         double mouseX, double mouseY) {
    const double dx = second.x - first.x;
    const double dy = second.y - first.y;
    const double lengthSquared = dx * dx + dy * dy;
    const double fraction = lengthSquared > 0.0
        ? std::clamp(((mouseX - first.x) * dx + (mouseY - first.y) * dy)
            / lengthSquared, 0.0, 1.0)
        : 0.0;
    Plot3DProbe result;
    result.traceIndex = traceIndex;
    result.segmentIndex = segmentIndex;
    result.segmentFraction = fraction;
    result.value = interpolate(firstValue, secondValue, fraction);
    result.projected = {
        first.x + fraction * dx,
        first.y + fraction * dy,
        first.depth + fraction * (second.depth - first.depth), true};
    result.distancePixels = std::hypot(result.projected.x - mouseX,
                                       result.projected.y - mouseY);
    return result;
}

}  // namespace

bool Plot3DCameraState::operator==(const Plot3DCameraState& other) const {
    return azimuthDegrees == other.azimuthDegrees
        && elevationDegrees == other.elevationDegrees
        && zoom == other.zoom && panX == other.panX && panY == other.panY
        && projection == other.projection;
}

bool Plot3DControlBounds::contains(double x, double y) const {
    return x >= left && x <= left + width && y >= top && y <= top + height;
}

Plot3DCamera::Plot3DCamera(Plot3DCameraState state) : state_(state) {
    normalize();
}

void Plot3DCamera::setState(Plot3DCameraState state) {
    state_ = state;
    normalize();
}

void Plot3DCamera::reset() { state_ = {}; }

void Plot3DCamera::normalize() {
    state_.azimuthDegrees = finiteOr(state_.azimuthDegrees, 42.0);
    state_.azimuthDegrees = std::fmod(state_.azimuthDegrees, 360.0);
    if (state_.azimuthDegrees < 0.0) state_.azimuthDegrees += 360.0;
    state_.elevationDegrees = std::clamp(
        finiteOr(state_.elevationDegrees, 24.0), -89.0, 89.0);
    state_.zoom = std::clamp(finiteOr(state_.zoom, 1.0), 0.12, 16.0);
    state_.panX = std::clamp(finiteOr(state_.panX, 0.0), -4.0, 4.0);
    state_.panY = std::clamp(finiteOr(state_.panY, 0.0), -4.0, 4.0);
}

void Plot3DCamera::orbit(double pixelDeltaX, double pixelDeltaY,
                         bool constrainAxis) {
    if (constrainAxis) {
        if (std::abs(pixelDeltaX) >= std::abs(pixelDeltaY)) pixelDeltaY = 0.0;
        else pixelDeltaX = 0.0;
    }
    state_.azimuthDegrees += pixelDeltaX * 0.35;
    state_.elevationDegrees -= pixelDeltaY * 0.35;
    normalize();
}

void Plot3DCamera::pan(double pixelDeltaX, double pixelDeltaY,
                       const Plot3DViewport& viewport) {
    const double side = std::max(1.0, std::min(viewport.width, viewport.height));
    state_.panX += pixelDeltaX / (side * 0.42 * state_.zoom);
    state_.panY -= pixelDeltaY / (side * 0.42 * state_.zoom);
    normalize();
}

void Plot3DCamera::zoomBy(double factor) {
    if (std::isfinite(factor) && factor > 0.0) state_.zoom *= factor;
    normalize();
}

void Plot3DCamera::zoomDrag(double pixelDeltaY) {
    zoomBy(std::exp(-pixelDeltaY * 0.012));
}

void Plot3DCamera::toggleProjection() {
    state_.projection = state_.projection == Plot3DProjection::perspective
        ? Plot3DProjection::orthographic : Plot3DProjection::perspective;
}

void Plot3DCamera::apply(Plot3DControl control) {
    switch (control) {
    case Plot3DControl::orbit_left: orbit(-22.0, 0.0); break;
    case Plot3DControl::orbit_right: orbit(22.0, 0.0); break;
    case Plot3DControl::orbit_up: orbit(0.0, -22.0); break;
    case Plot3DControl::orbit_down: orbit(0.0, 22.0); break;
    case Plot3DControl::pan_left: state_.panX -= 0.09; break;
    case Plot3DControl::pan_right: state_.panX += 0.09; break;
    case Plot3DControl::pan_up: state_.panY += 0.09; break;
    case Plot3DControl::pan_down: state_.panY -= 0.09; break;
    case Plot3DControl::zoom_in: zoomBy(1.18); break;
    case Plot3DControl::zoom_out: zoomBy(1.0 / 1.18); break;
    case Plot3DControl::toggle_projection: toggleProjection(); break;
    case Plot3DControl::reset: reset(); break;
    case Plot3DControl::none: break;
    }
    normalize();
}

Plot3DProjectedPoint Plot3DCamera::project(
    Plot3DPoint point, const Plot3DBounds& bounds,
    const Plot3DViewport& viewport) const {
    const Plot3DPoint center{
        (bounds.minimum.x + bounds.maximum.x) * 0.5,
        (bounds.minimum.y + bounds.maximum.y) * 0.5,
        (bounds.minimum.z + bounds.maximum.z) * 0.5};
    const double span = std::max({
        bounds.maximum.x - bounds.minimum.x,
        bounds.maximum.y - bounds.minimum.y,
        bounds.maximum.z - bounds.minimum.z, 1.0e-14});
    Plot3DPoint normalized = subtract(point, center);
    normalized.x *= 2.0 / span;
    normalized.y *= 2.0 / span;
    normalized.z *= 2.0 / span;

    const double azimuth = radians(state_.azimuthDegrees);
    const double elevation = radians(state_.elevationDegrees);
    const Plot3DPoint right{-std::sin(azimuth), std::cos(azimuth), 0.0};
    const Plot3DPoint up{-std::sin(elevation) * std::cos(azimuth),
                         -std::sin(elevation) * std::sin(azimuth),
                         std::cos(elevation)};
    const Plot3DPoint forward{std::cos(elevation) * std::cos(azimuth),
                              std::cos(elevation) * std::sin(azimuth),
                              std::sin(elevation)};
    const double horizontal = dot(normalized, right);
    const double vertical = dot(normalized, up);
    const double depth = dot(normalized, forward);
    double perspective = 1.0;
    if (state_.projection == Plot3DProjection::perspective)
        perspective = 3.4 / std::max(0.45, 3.4 - depth);
    const double side = std::max(1.0, std::min(viewport.width, viewport.height));
    const double scale = side * 0.39 * state_.zoom;
    return {viewport.left + viewport.width * 0.5
                + (horizontal * perspective + state_.panX) * scale,
            viewport.top + viewport.height * 0.5
                - (vertical * perspective + state_.panY) * scale,
            depth, std::isfinite(horizontal) && std::isfinite(vertical)};
}

Plot3DBounds plot3DBounds(const Figure& figure) {
    Plot3DBounds result{{std::numeric_limits<double>::infinity(),
                         std::numeric_limits<double>::infinity(),
                         std::numeric_limits<double>::infinity()},
                        {-std::numeric_limits<double>::infinity(),
                         -std::numeric_limits<double>::infinity(),
                         -std::numeric_limits<double>::infinity()}};
    for (const auto& trace : figure.traces) {
        const auto count = std::min({trace.x.size(), trace.y.size(), trace.z.size()});
        for (std::size_t index = 0; index < count; ++index) {
            const auto point = tracePoint(trace, index);
            if (!finite(point)) continue;
            result.minimum.x = std::min(result.minimum.x, point.x);
            result.minimum.y = std::min(result.minimum.y, point.y);
            result.minimum.z = std::min(result.minimum.z, point.z);
            result.maximum.x = std::max(result.maximum.x, point.x);
            result.maximum.y = std::max(result.maximum.y, point.y);
            result.maximum.z = std::max(result.maximum.z, point.z);
        }
    }
    if (!finite(result.minimum) || !finite(result.maximum)) return {};
    auto expand = [](double& low, double& high) {
        if (std::abs(high - low) < 1.0e-14) { low -= 0.5; high += 0.5; }
        const double margin = (high - low) * 0.045;
        low -= margin; high += margin;
    };
    expand(result.minimum.x, result.maximum.x);
    expand(result.minimum.y, result.maximum.y);
    expand(result.minimum.z, result.maximum.z);
    return result;
}

std::optional<Plot3DProbe> nearestPlot3DProbe(
    const Figure& figure, const Plot3DCamera& camera,
    const Plot3DViewport& viewport, double pixelX, double pixelY,
    double maximumDistancePixels, std::optional<std::size_t> onlyTrace) {
    const auto bounds = plot3DBounds(figure);
    std::optional<Plot3DProbe> best;
    for (std::size_t traceIndex = 0; traceIndex < figure.traces.size(); ++traceIndex) {
        if (onlyTrace && *onlyTrace != traceIndex) continue;
        const auto& trace = figure.traces[traceIndex];
        const auto count = std::min({trace.x.size(), trace.y.size(), trace.z.size()});
        if (!count) continue;
        Plot3DPoint previousValue = tracePoint(trace, 0);
        auto previous = camera.project(previousValue, bounds, viewport);
        for (std::size_t index = 1; index < count; ++index) {
            const auto currentValue = tracePoint(trace, index);
            const auto current = camera.project(currentValue, bounds, viewport);
            if (finite(previousValue) && finite(currentValue)
                && previous.visible && current.visible) {
                const auto candidate = segmentProbe(previousValue, currentValue,
                    previous, current, traceIndex, index - 1, pixelX, pixelY);
                if (!best || candidate.distancePixels < best->distancePixels)
                    best = candidate;
            }
            previousValue = currentValue;
            previous = current;
        }
    }
    if (!best || best->distancePixels > maximumDistancePixels) return std::nullopt;
    return best;
}

std::optional<Plot3DProbe> resolvePlot3DProbe(
    const Figure& figure, const Plot3DCamera& camera,
    const Plot3DViewport& viewport, std::size_t traceIndex,
    std::size_t segmentIndex, double segmentFraction) {
    if (traceIndex >= figure.traces.size()) return std::nullopt;
    const auto& trace = figure.traces[traceIndex];
    const auto count = std::min({trace.x.size(), trace.y.size(), trace.z.size()});
    if (!count) return std::nullopt;
    const auto firstIndex = std::min(segmentIndex, count - 1);
    const auto secondIndex = std::min(firstIndex + 1, count - 1);
    const auto firstValue = tracePoint(trace, firstIndex);
    const auto secondValue = tracePoint(trace, secondIndex);
    const auto bounds = plot3DBounds(figure);
    if (!finite(firstValue) || !finite(secondValue)) return std::nullopt;
    Plot3DProbe result;
    result.traceIndex = traceIndex;
    result.segmentIndex = firstIndex;
    result.segmentFraction = std::clamp(segmentFraction, 0.0, 1.0);
    result.value = interpolate(firstValue, secondValue, result.segmentFraction);
    result.projected = camera.project(result.value, bounds, viewport);
    result.distancePixels = 0.0;
    return result.projected.visible ? std::optional<Plot3DProbe>(result)
                                    : std::nullopt;
}

std::vector<Plot3DControlBounds> plot3DControls(const Plot3DViewport& viewport) {
    constexpr double size = 40.0;
    constexpr double gap = 4.0;
    const double total = size * 4.0 + gap * 3.0;
    const double left = std::max(viewport.left + 6.0,
        viewport.left + viewport.width - total - 8.0);
    const double top = viewport.top + 8.0;
    const Plot3DControl rows[3][4] = {
        {Plot3DControl::orbit_left, Plot3DControl::orbit_up,
         Plot3DControl::orbit_down, Plot3DControl::orbit_right},
        {Plot3DControl::pan_left, Plot3DControl::pan_up,
         Plot3DControl::pan_down, Plot3DControl::pan_right},
        {Plot3DControl::zoom_in, Plot3DControl::zoom_out,
         Plot3DControl::toggle_projection, Plot3DControl::reset}};
    std::vector<Plot3DControlBounds> result;
    result.reserve(12);
    for (int row = 0; row < 3; ++row)
        for (int column = 0; column < 4; ++column)
            result.push_back({rows[row][column], left + column * (size + gap),
                              top + row * (size + gap), size, size});
    return result;
}

Plot3DControl hitPlot3DControl(const Plot3DViewport& viewport,
                               double x, double y) {
    for (const auto& control : plot3DControls(viewport))
        if (control.contains(x, y)) return control.control;
    return Plot3DControl::none;
}

std::string plot3DControlName(Plot3DControl control) {
    switch (control) {
    case Plot3DControl::orbit_left: return "Rotate left";
    case Plot3DControl::orbit_right: return "Rotate right";
    case Plot3DControl::orbit_up: return "Rotate up";
    case Plot3DControl::orbit_down: return "Rotate down";
    case Plot3DControl::pan_left: return "Move left";
    case Plot3DControl::pan_right: return "Move right";
    case Plot3DControl::pan_up: return "Move up";
    case Plot3DControl::pan_down: return "Move down";
    case Plot3DControl::zoom_in: return "Zoom in";
    case Plot3DControl::zoom_out: return "Zoom out";
    case Plot3DControl::toggle_projection: return "Toggle perspective";
    case Plot3DControl::reset: return "Reset view";
    case Plot3DControl::none: return {};
    }
    return {};
}

}  // namespace dsl
