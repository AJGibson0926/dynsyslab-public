#include "dynsyslab_viewer/Plot3D.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

bool near(double first, double second, double tolerance = 1.0e-8) {
    return std::abs(first - second) <= tolerance;
}

}  // namespace

int main() {
    try {
        const dsl::Plot3DViewport viewport{0, 0, 800, 600};
        const dsl::Plot3DBounds bounds{{-1, -2, -3}, {1, 2, 3}};
        dsl::Plot3DCamera camera;
        const auto initial = camera.state();
        const auto first = camera.project({1, 0, 0}, bounds, viewport);
        require(first.visible && std::isfinite(first.x) && std::isfinite(first.y),
                "default camera did not produce a finite projection");

        camera.orbit(60, -30);
        const auto rotated = camera.project({1, 0, 0}, bounds, viewport);
        require(!near(first.x, rotated.x) || !near(first.y, rotated.y),
                "orbit did not rotate the projected point");
        camera.pan(40, -20, viewport);
        require(camera.state().panX > 0 && camera.state().panY > 0,
                "pan did not move the camera frame");
        camera.zoomBy(2.0);
        require(camera.state().zoom > initial.zoom,
                "zoom did not change camera magnification");
        camera.toggleProjection();
        for (int i=0;i<2;++i) {
            const auto anchor=camera.project({0.3,-0.2,0.1},bounds,viewport);
            camera.zoomAt(1.12, anchor.x, anchor.y, viewport);
            const auto after=camera.project({0.3,-0.2,0.1},bounds,viewport);
            require(near(after.x,anchor.x) && near(after.y,anchor.y),"wheel shifted the 3D point under the cursor");
            camera.toggleProjection();
        }
        require(camera.state().projection == dsl::Plot3DProjection::orthographic,
                "projection toggle did not select orthographic mode");
        camera.orbit(0, -10000);
        camera.zoomBy(1.0e10);
        require(camera.state().elevationDegrees <= 89.0
            && camera.state().zoom <= 16.0,
            "camera limits did not clamp unstable input");
        camera.reset();
        require(camera.state() == initial, "camera reset was not deterministic");

        dsl::Figure figure;
        figure.kind = "three_d";
        dsl::Trace trace;
        trace.x = {-1, 0, 1}; trace.y = {-1, 0, 1}; trace.z = {-1, 0, 1};
        figure.traces.push_back(trace);
        const auto dataBounds = dsl::plot3DBounds(figure);
        require(dataBounds.minimum.x < -1 && dataBounds.maximum.z > 1,
                "three-dimensional bounds omitted their display margin");
        const auto midpoint = camera.project({0, 0, 0}, dataBounds, viewport);
        const auto probe = dsl::nearestPlot3DProbe(figure, camera, viewport,
            midpoint.x, midpoint.y, 2.0);
        require(probe && near(probe->value.x, 0.0, 1.0e-6)
            && near(probe->value.z, 0.0, 1.0e-6),
            "3D curve probe did not preserve the source coordinates");
        const auto restored = dsl::resolvePlot3DProbe(figure, camera, viewport,
            probe->traceIndex, probe->segmentIndex, probe->segmentFraction);
        require(restored && near(restored->projected.x, probe->projected.x, 1.0e-6),
                "3D curve probe could not be restored in the camera frame");

        const auto controls = dsl::plot3DControls(viewport);
        require(controls.size() == 12, "3D gizmo layout is incomplete");
        for (const auto& control : controls)
            require(control.width >= 40 && control.height >= 40,
                    "3D gizmo target is smaller than the design-system target");
        require(dsl::hitPlot3DControl(viewport,
            controls.front().left + 10, controls.front().top + 10)
                == dsl::Plot3DControl::orbit_left,
                "3D gizmo hit testing failed");

        std::cout << "dynsyslab_plot_3d_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_plot_3d_tests: FAIL: " << error.what() << '\n';
        return 1;
    }
}
