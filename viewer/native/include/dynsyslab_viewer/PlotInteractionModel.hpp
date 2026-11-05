#pragma once
#include "dynsyslab_viewer/PlotScale.hpp"

#include "dynsyslab_viewer/Model.hpp"

#include <cstddef>
#include <optional>

namespace dsl {

struct PlotTransform {
    double left{};
    double top{};
    double width{1.0};
    double height{1.0};
    double xMin{-1.0};
    double xMax{1.0};
    double yMin{-1.0};
    double yMax{1.0};
    ValueScale yScale;
};

double plotPixelY(const PlotTransform& transform, double value);
double plotDataY(const PlotTransform& transform, double pixel);

struct CurveProbe {
    std::size_t traceIndex{};
    std::size_t segmentIndex{};
    double segmentFraction{};
    double x{};
    double y{};
    double pixelX{};
    double pixelY{};
    double distancePixels{};
    std::optional<double> z;
};

enum class PlotPointerKind { background, curve, pinned_probe };
enum class PlotGestureKind { none, probe_drag, background_pan };

struct PlotPointerDecision {
    PlotPointerKind kind{PlotPointerKind::background};
    std::optional<CurveProbe> probe;
};

struct PlotInteractionState {
    PlotTransform transform{};
    std::optional<CurveProbe> hover;
    std::optional<CurveProbe> pinned;
    PlotPointerKind pointer{PlotPointerKind::background};
    PlotGestureKind gesture{PlotGestureKind::none};
};

[[nodiscard]] bool containsPixel(const PlotTransform& transform,
                                 double pixelX, double pixelY);
[[nodiscard]] std::optional<CurveProbe> nearestCurveProbe(
    const Figure& figure, const PlotTransform& transform,
    double pixelX, double pixelY, double maximumDistancePixels,
    std::optional<std::size_t> traceIndex = std::nullopt);
[[nodiscard]] std::optional<CurveProbe> resolveCurveProbe(
    const Figure& figure, const PlotTransform& transform,
    std::size_t traceIndex, std::size_t segmentIndex, double segmentFraction);
[[nodiscard]] PlotTransform panPlotTransform(const PlotTransform& transform,
                                             double pixelDeltaX,
                                             double pixelDeltaY);
[[nodiscard]] PlotPointerDecision inspectPlotPointer(
    const Figure& figure, const PlotTransform& transform,
    const std::optional<CurveProbe>& pinned, double pixelX, double pixelY,
    double hitRadiusPixels = 9.0);

class PlotInteractionSession {
public:
    PlotInteractionSession(const Figure& figure, PlotTransform transform);
    void move(double pixelX, double pixelY);
    void press(double pixelX, double pixelY);
    void drag(double pixelX, double pixelY);
    void release();
    [[nodiscard]] const PlotInteractionState& state() const { return state_; }

private:
    const Figure* figure_{};
    PlotInteractionState state_{};
    PlotTransform gestureStart_{};
    double pressX_{};
    double pressY_{};
};

}  // namespace dsl
