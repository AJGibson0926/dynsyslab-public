#include "dynsyslab_viewer/PlotInteractionModel.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace dsl {
double plotPixelY(const PlotTransform& t, double value) {
    const auto low = t.yScale.forward(t.yMin), high = t.yScale.forward(t.yMax);
    return t.top+t.height-(t.yScale.forward(value)-low)/(high-low)*t.height;
}
double plotDataY(const PlotTransform& t, double pixel) {
    const auto low = t.yScale.forward(t.yMin), high = t.yScale.forward(t.yMax);
    return t.yScale.inverse(low+(t.top+t.height-pixel)/t.height*(high-low));
}
namespace {

struct Sample {
    double x{};
    double y{};
    double pixelX{};
    double pixelY{};
    bool valid{};
};

Sample sample(const Trace& trace, std::size_t index,
              const PlotTransform& transform) {
    if (index >= trace.x.size() || index >= trace.y.size()) return {};
    double x = trace.x[index], y = trace.y[index];
    if (!trace.z.empty() && index < trace.z.size()) {
        x -= 0.32 * trace.z[index];
        y += 0.20 * trace.z[index];
    }
    const double xRange = transform.xMax - transform.xMin;
    const double yRange = transform.yMax - transform.yMin;
    if (!std::isfinite(x) || !std::isfinite(y) || xRange == 0.0 || yRange == 0.0)
        return {};
    return {x, y,
        transform.left + (x - transform.xMin) / xRange * transform.width,
        plotPixelY(transform,y),
        true};
}

CurveProbe pointProbe(const Sample& value, std::size_t traceIndex,
                      std::size_t segmentIndex, double fraction,
                      double mouseX, double mouseY) {
    const double dx = value.pixelX - mouseX;
    const double dy = value.pixelY - mouseY;
    return {traceIndex, segmentIndex, fraction, value.x, value.y,
            value.pixelX, value.pixelY, std::hypot(dx, dy)};
}

CurveProbe segmentProbe(const Sample& first, const Sample& second,
                        std::size_t traceIndex, std::size_t segmentIndex,
                        double mouseX, double mouseY, const PlotTransform& transform) {
    const double dx = second.pixelX - first.pixelX;
    const double dy = second.pixelY - first.pixelY;
    const double lengthSquared = dx * dx + dy * dy;
    const double fraction = lengthSquared > 0.0
        ? std::clamp(((mouseX - first.pixelX) * dx
            + (mouseY - first.pixelY) * dy) / lengthSquared, 0.0, 1.0)
        : 0.0;
    Sample projected;
    projected.x = first.x + fraction * (second.x - first.x);
    projected.y = transform.yScale.inverse(transform.yScale.forward(first.y)
        + fraction * (transform.yScale.forward(second.y)-transform.yScale.forward(first.y)));
    projected.pixelX = first.pixelX + fraction * dx;
    projected.pixelY = first.pixelY + fraction * dy;
    projected.valid = true;
    return pointProbe(projected, traceIndex, segmentIndex, fraction, mouseX, mouseY);
}

}  // namespace

bool containsPixel(const PlotTransform& transform, double pixelX, double pixelY) {
    return pixelX >= transform.left && pixelX <= transform.left + transform.width
        && pixelY >= transform.top && pixelY <= transform.top + transform.height;
}

std::optional<CurveProbe> nearestCurveProbe(
    const Figure& figure, const PlotTransform& transform,
    double pixelX, double pixelY, double maximumDistancePixels,
    std::optional<std::size_t> requestedTrace) {
    std::optional<CurveProbe> best;
    for (std::size_t traceIndex = 0; traceIndex < figure.traces.size(); ++traceIndex) {
        if (requestedTrace && *requestedTrace != traceIndex) continue;
        const Trace& trace = figure.traces[traceIndex];
        if (trace.heatmap) continue;
        const std::size_t count = std::min(trace.x.size(), trace.y.size());
        if (!count) continue;
        Sample previous = sample(trace, 0, transform);
        if (previous.valid) {
            const auto candidate = pointProbe(previous, traceIndex, 0, 0.0, pixelX, pixelY);
            if (!best || candidate.distancePixels < best->distancePixels) best = candidate;
        }
        for (std::size_t index = 1; index < count; ++index) {
            const Sample current = sample(trace, index, transform);
            if (previous.valid && current.valid) {
                const auto candidate = segmentProbe(previous, current, traceIndex,
                                                    index - 1, pixelX, pixelY, transform);
                if (!best || candidate.distancePixels < best->distancePixels) best = candidate;
            }
            previous = current;
        }
    }
    if (!best || best->distancePixels > maximumDistancePixels) return std::nullopt;
    return best;
}

std::optional<CurveProbe> resolveCurveProbe(
    const Figure& figure, const PlotTransform& transform,
    std::size_t traceIndex, std::size_t segmentIndex, double segmentFraction) {
    if (traceIndex >= figure.traces.size()) return std::nullopt;
    const Trace& trace = figure.traces[traceIndex];
    const std::size_t count = std::min(trace.x.size(), trace.y.size());
    if (!count) return std::nullopt;
    const std::size_t firstIndex = std::min(segmentIndex, count - 1);
    const Sample first = sample(trace, firstIndex, transform);
    if (!first.valid) return std::nullopt;
    if (firstIndex + 1 >= count)
        return pointProbe(first, traceIndex, firstIndex, 0.0,
                          first.pixelX, first.pixelY);
    const Sample second = sample(trace, firstIndex + 1, transform);
    if (!second.valid) return std::nullopt;
    const double fraction = std::clamp(segmentFraction, 0.0, 1.0);
    Sample projected;
    projected.x = first.x + fraction * (second.x - first.x);
    projected.y = transform.yScale.inverse(transform.yScale.forward(first.y)
        + fraction * (transform.yScale.forward(second.y)-transform.yScale.forward(first.y)));
    projected.pixelX = first.pixelX + fraction * (second.pixelX - first.pixelX);
    projected.pixelY = first.pixelY + fraction * (second.pixelY - first.pixelY);
    projected.valid = true;
    return pointProbe(projected, traceIndex, firstIndex, fraction,
                      projected.pixelX, projected.pixelY);
}

std::optional<AxisLimits> zoomPlotLimits(const PlotTransform& t, double pixelX,
                                         double pixelY, double rangeFactor) {
    if (!(t.width > 0 && t.height > 0 && rangeFactor > 0) || !std::isfinite(rangeFactor)) return {};
    const double fx = std::clamp((pixelX - t.left) / t.width, 0.0, 1.0);
    const double fy = std::clamp((t.top + t.height - pixelY) / t.height, 0.0, 1.0);
    const double yLow = t.yScale.forward(t.yMin), yHigh = t.yScale.forward(t.yMax);
    const double dx = (1.0 - rangeFactor) * (t.xMax - t.xMin);
    const double dy = (1.0 - rangeFactor) * (yHigh - yLow);
    AxisLimits limits{t.xMin + fx * dx, t.xMax - (1.0 - fx) * dx,
        t.yScale.inverse(yLow + fy * dy), t.yScale.inverse(yHigh - (1.0 - fy) * dy)};
    if (!(limits.xMin < limits.xMax && limits.yMin < limits.yMax)
        || !std::isfinite(limits.xMax - limits.xMin)
        || !std::isfinite(t.yScale.forward(limits.yMax) - t.yScale.forward(limits.yMin))) return {};
    return limits;
}

PlotTransform panPlotTransform(const PlotTransform& transform,
                               double pixelDeltaX, double pixelDeltaY) {
    PlotTransform result = transform;
    const double xShift = -pixelDeltaX / std::max(transform.width, 1.0)
        * (transform.xMax - transform.xMin);
    const double yShift = pixelDeltaY / std::max(transform.height, 1.0)
        * (transform.yScale.forward(transform.yMax) - transform.yScale.forward(transform.yMin));
    result.xMin += xShift; result.xMax += xShift;
    result.yMin = transform.yScale.inverse(transform.yScale.forward(transform.yMin)+yShift);
    result.yMax = transform.yScale.inverse(transform.yScale.forward(transform.yMax)+yShift);
    return result;
}

PlotPointerDecision inspectPlotPointer(
    const Figure& figure, const PlotTransform& transform,
    const std::optional<CurveProbe>& pinned, double pixelX, double pixelY,
    double hitRadiusPixels) {
    if (pinned && std::hypot(pinned->pixelX - pixelX,
                             pinned->pixelY - pixelY) <= hitRadiusPixels)
        return {PlotPointerKind::pinned_probe, pinned};
    const auto curve = nearestCurveProbe(figure, transform, pixelX, pixelY,
                                         hitRadiusPixels);
    return curve ? PlotPointerDecision{PlotPointerKind::curve, curve}
                 : PlotPointerDecision{};
}

PlotInteractionSession::PlotInteractionSession(
    const Figure& figure, PlotTransform transform) : figure_(&figure) {
    state_.transform = transform;
}

void PlotInteractionSession::move(double pixelX, double pixelY) {
    if (state_.gesture != PlotGestureKind::none) {
        drag(pixelX, pixelY); return;
    }
    const auto decision = inspectPlotPointer(*figure_, state_.transform,
        state_.pinned, pixelX, pixelY);
    state_.pointer = decision.kind;
    state_.hover = decision.kind == PlotPointerKind::curve
        ? decision.probe : std::nullopt;
}

void PlotInteractionSession::press(double pixelX, double pixelY) {
    const auto decision = inspectPlotPointer(*figure_, state_.transform,
        state_.pinned, pixelX, pixelY);
    if (decision.kind == PlotPointerKind::curve && decision.probe)
        state_.pinned = decision.probe;
    if (decision.kind != PlotPointerKind::background) {
        state_.hover = state_.pinned;
        state_.gesture = PlotGestureKind::probe_drag;
        state_.pointer = PlotPointerKind::pinned_probe;
        return;
    }
    state_.gesture = PlotGestureKind::background_pan;
    state_.pointer = PlotPointerKind::background;
    gestureStart_ = state_.transform; pressX_ = pixelX; pressY_ = pixelY;
    state_.hover.reset();
}

void PlotInteractionSession::drag(double pixelX, double pixelY) {
    if (state_.gesture == PlotGestureKind::probe_drag && state_.pinned) {
        const auto probe = nearestCurveProbe(*figure_, state_.transform,
            pixelX, pixelY, std::numeric_limits<double>::infinity(),
            state_.pinned->traceIndex);
        if (probe) state_.pinned = probe;
        state_.hover = state_.pinned;
    } else if (state_.gesture == PlotGestureKind::background_pan) {
        state_.transform = panPlotTransform(gestureStart_,
            pixelX - pressX_, pixelY - pressY_);
        if (state_.pinned) state_.pinned = resolveCurveProbe(*figure_,
            state_.transform, state_.pinned->traceIndex,
            state_.pinned->segmentIndex, state_.pinned->segmentFraction);
        state_.hover.reset();
    }
}

void PlotInteractionSession::release() {
    state_.gesture = PlotGestureKind::none;
}

}  // namespace dsl
