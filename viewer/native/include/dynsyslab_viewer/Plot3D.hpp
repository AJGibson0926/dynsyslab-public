#pragma once

#include "dynsyslab_viewer/Model.hpp"

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace dsl {

struct Plot3DPoint {
    double x{};
    double y{};
    double z{};
};

struct Plot3DBounds {
    Plot3DPoint minimum{-1.0, -1.0, -1.0};
    Plot3DPoint maximum{1.0, 1.0, 1.0};
};

struct Plot3DViewport {
    double left{};
    double top{};
    double width{1.0};
    double height{1.0};
};

enum class Plot3DProjection { perspective, orthographic };

struct Plot3DCameraState {
    double azimuthDegrees = 42.0;
    double elevationDegrees = 24.0;
    double zoom = 1.0;
    double panX = 0.0;
    double panY = 0.0;
    Plot3DProjection projection = Plot3DProjection::perspective;
    bool boxAxes = false;  // independent per-axis scaling (stem plots)

    bool operator==(const Plot3DCameraState& other) const;
    bool operator!=(const Plot3DCameraState& other) const { return !(*this == other); }

    // Camera a figure opens with (and resets to); defaults unless it sets one.
    template <class FigureT>
    static Plot3DCameraState initialFor(const FigureT& figure) {
        Plot3DCameraState state;
        if (figure.cameraAzimuthDegrees) state.azimuthDegrees = *figure.cameraAzimuthDegrees;
        if (figure.cameraElevationDegrees) state.elevationDegrees = *figure.cameraElevationDegrees;
        state.boxAxes = figure.boxAxes3D;
        if (state.boxAxes) {  // a full cube must fit the panel; no perspective distortion of heights
            state.zoom = 0.7;
            state.projection = Plot3DProjection::orthographic;
        }
        return state;
    }
};

struct Plot3DProjectedPoint {
    double x{};
    double y{};
    double depth{};
    bool visible{};
};

struct Plot3DProbe {
    std::size_t traceIndex{};
    std::size_t segmentIndex{};
    double segmentFraction{};
    Plot3DPoint value{};
    Plot3DProjectedPoint projected{};
    double distancePixels{};
};

enum class Plot3DControl {
    none,
    orbit_left,
    orbit_up,
    orbit_down,
    orbit_right,
    pan_left,
    pan_up,
    pan_down,
    pan_right,
    zoom_in,
    zoom_out,
    toggle_projection,
    reset
};

struct Plot3DControlBounds {
    Plot3DControl control{Plot3DControl::none};
    double left{};
    double top{};
    double width{};
    double height{};

    [[nodiscard]] bool contains(double x, double y) const;
};

class Plot3DCamera {
public:
    Plot3DCamera() = default;
    explicit Plot3DCamera(Plot3DCameraState state);

    [[nodiscard]] const Plot3DCameraState& state() const { return state_; }
    void setState(Plot3DCameraState state);
    void reset();
    void orbit(double pixelDeltaX, double pixelDeltaY, bool constrainAxis = false);
    void pan(double pixelDeltaX, double pixelDeltaY, const Plot3DViewport& viewport);
    void zoomBy(double factor);
    void zoomAt(double factor, double pixelX, double pixelY, const Plot3DViewport& viewport);
    void zoomDrag(double pixelDeltaY);
    void toggleProjection();
    void apply(Plot3DControl control);

    [[nodiscard]] Plot3DProjectedPoint project(
        Plot3DPoint point, const Plot3DBounds& bounds,
        const Plot3DViewport& viewport) const;

private:
    Plot3DCameraState state_{};
    void normalize();
};

// Height (z) scale of box-axis 3-D figures (stem plots), chosen by the figure's
// valueScale: "linear", "log"/"symlog" (logarithmic heights), or "auto" (log when the
// positive heights span at least three decades). Stems start at zero, so on a log
// axis every height at or below the floor (at most six decades under the tallest)
// is drawn at the floor. Trajectory 3-D figures always stay linear.
struct Plot3DHeightScale {
    bool logarithmic = false;
    double floorExponent = 0.0;          // log10 of the floor value
    [[nodiscard]] double forward(double z) const;
    [[nodiscard]] double inverse(double z) const;
};
[[nodiscard]] Plot3DHeightScale plot3DHeightScale(const Figure& figure);
// Copy of the figure with z mapped through its height scale (identity when linear).
[[nodiscard]] Figure plot3DDisplayFigure(const Figure& figure);

[[nodiscard]] Plot3DBounds plot3DBounds(const Figure& figure);
[[nodiscard]] std::optional<Plot3DProbe> nearestPlot3DProbe(
    const Figure& figure, const Plot3DCamera& camera,
    const Plot3DViewport& viewport, double pixelX, double pixelY,
    double maximumDistancePixels,
    std::optional<std::size_t> traceIndex = std::nullopt);
[[nodiscard]] std::optional<Plot3DProbe> resolvePlot3DProbe(
    const Figure& figure, const Plot3DCamera& camera,
    const Plot3DViewport& viewport, std::size_t traceIndex,
    std::size_t segmentIndex, double segmentFraction);

[[nodiscard]] std::vector<Plot3DControlBounds> plot3DControls(
    const Plot3DViewport& viewport);
[[nodiscard]] Plot3DControl hitPlot3DControl(
    const Plot3DViewport& viewport, double x, double y);
[[nodiscard]] std::string plot3DControlName(Plot3DControl control);

}  // namespace dsl
