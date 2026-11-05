#include "dynsyslab_viewer/MathText.hpp"
#include "dynsyslab_viewer/PlotInteractionModel.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

bool near(double first, double second) {
    return std::abs(first - second) < 1.0e-9;
}

}  // namespace

int main() {
    try {
        dsl::Figure figure;
        dsl::Trace curve;
        curve.x = {0.0, 10.0}; curve.y = {0.0, 10.0};
        figure.traces.push_back(curve);
        dsl::Trace heat; heat.heatmap = true;
        heat.x = {0.0, 10.0}; heat.y = {0.0, 10.0};
        figure.traces.push_back(heat);
        const dsl::PlotTransform transform{0.0, 0.0, 100.0, 100.0,
                                           0.0, 10.0, 0.0, 10.0};

        require(dsl::containsPixel(transform, 25.0, 75.0),
            "plot-background containment failed");
        const auto probe = dsl::nearestCurveProbe(figure, transform,
                                                   50.0, 50.0, 2.0);
        require(probe.has_value(), "curve segment was not hit");
        require(probe->traceIndex == 0 && near(probe->x, 5.0)
            && near(probe->y, 5.0), "curve hit did not interpolate data coordinates");
        require(!dsl::nearestCurveProbe(figure, transform, 50.0, 75.0, 3.0),
            "distant background incorrectly hit a curve");
        const auto restored = dsl::resolveCurveProbe(figure, transform,
            probe->traceIndex, probe->segmentIndex, probe->segmentFraction);
        require(restored && near(restored->pixelX, 50.0),
            "pinned probe could not be restored after rendering");

        const auto panned = dsl::panPlotTransform(transform, 10.0, -20.0);
        require(near(panned.xMin, -1.0) && near(panned.xMax, 9.0)
            && near(panned.yMin, -2.0) && near(panned.yMax, 8.0),
            "background drag did not pan both axes correctly");

        dsl::PlotInteractionSession session(figure, transform);
        session.move(50.0, 50.0);
        require(session.state().pointer == dsl::PlotPointerKind::curve
            && session.state().hover, "session did not enter curve-pointer mode");
        session.press(50.0, 50.0);
        require(session.state().pinned
            && session.state().gesture == dsl::PlotGestureKind::probe_drag,
            "session did not pin the curve probe");
        session.drag(75.0, 25.0); session.release();
        require(session.state().pinned && near(session.state().pinned->x, 7.5),
            "session did not drag the pinned probe along its curve");
        session.move(10.0, 20.0);
        require(session.state().pointer == dsl::PlotPointerKind::background,
            "session did not enter background-move mode");
        session.press(10.0, 20.0); session.drag(20.0, 30.0); session.release();
        require(near(session.state().transform.xMin, -1.0)
            && near(session.state().transform.yMin, 1.0),
            "session background drag did not pan the plot");

        const auto spans = dsl::parseMathText(L"$\\psi_{12}^{*}$");
        require(spans.size() == 3, "legend math-text span count failed");
        require(spans[0].text == L"\u03c8", "legend Greek symbol parsing failed");
        require(spans[1].text == L"12"
            && spans[1].role == dsl::MathTextRole::subscript,
            "legend subscript parsing failed");
        require(spans[2].text == L"*"
            && spans[2].role == dsl::MathTextRole::superscript,
            "legend superscript parsing failed");

        std::cout << "dynsyslab_plot_interaction_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_plot_interaction_tests: FAIL: "
                  << error.what() << '\n';
        return 1;
    }
}
