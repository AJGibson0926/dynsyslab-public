#include "PhaseContext.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
namespace dsl {
void appendPhaseSurfaces(Figure &figure, const json::Value &result, std::size_t n, double sx,
                         double sy) {
    if (figure.kind != "phase")
        return;
    for (const auto &mesh : result.at("manifold_surfaces").array()) {
        const auto status = mesh.at("status").string();
        if (status != "resolved" && status != "truncated")
            continue;
        const auto &vertices = mesh.at("vertices").array(),
                   &offsets = mesh.at("front_offsets").array();
        if (vertices.size() % n || offsets.size() < 3)
            continue;
        const auto count = vertices.size() / n;
        std::vector<std::size_t> starts;
        bool valid = true;
        for (const auto &q : offsets) {
            double v = q.number(-1);
            if (!std::isfinite(v) || v < 1 || v >= count || v != std::floor(v) ||
                (!starts.empty() && v <= starts.back())) {
                valid = false;
                break;
            }
            starts.push_back(std::size_t(v));
        }
        if (!valid)
            continue;
        starts.push_back(count);
        Trace trace;
        const bool stable = mesh.at("kind").string() == "stable_manifold_surface";
        trace.role = stable ? "phase_stable_surface" : "phase_unstable_surface";
        trace.label = std::string(stable ? "stable surface" : "unstable surface") +
                      (n > 2 ? " (projection)" : "");
        trace.color = stable ? "#9aaeb5" : "#b9a798";
        trace.lineStyle = "solid";
        auto point = [&](std::size_t i) {
            const double x = vertices[i * n].number(std::numeric_limits<double>::quiet_NaN()),
                         y = vertices[i * n + 1].number(std::numeric_limits<double>::quiet_NaN());
            trace.x.push_back(x * sx);
            trace.y.push_back(y * sy);
        };
        auto gap = [&]() {
            trace.x.push_back(std::numeric_limits<double>::quiet_NaN());
            trace.y.push_back(std::numeric_limits<double>::quiet_NaN());
        };
        // A few complete fronts convey surface geometry without opaque fills or dense triangles.
        const auto rings = starts.size() - 1;
        for (std::size_t k = 0; k < rings; k += std::max(std::size_t(1), (rings + 5) / 6)) {
            for (std::size_t i = starts[k]; i < starts[k + 1]; ++i)
                point(i);
            point(starts[k]);
            gap();
        }
        if (!trace.x.empty()) {
            if (figure.alternateCoordinates)
                figure.alternateCoordinates->traceLabels.push_back(trace.label);
            figure.traces.push_back(std::move(trace));
        }
    }
}
} // namespace dsl
