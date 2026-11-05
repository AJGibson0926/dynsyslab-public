#include "DisplayCoordinates.hpp"
#include "PhaseContext.hpp"
#include "dynsyslab_viewer/PlotPresentation.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
namespace dsl {
namespace {
bool finite(const json::Value &v) {
    return std::holds_alternative<double>(v.storage()) && std::isfinite(v.number());
}
void segment(Trace &trace, double x, double y, double u, double v) {
    const double gap = std::numeric_limits<double>::quiet_NaN();
    trace.x.insert(trace.x.end(), {x, u, gap});
    trace.y.insert(trace.y.end(), {y, v, gap});
}
} // namespace
void appendPhaseContext(Figure &figure, const json::Value &metadata, const json::Value &result) {
    if (figure.kind != "phase" || result.at("cache").at("status").string() == "miss")
        return;
    const bool periodic = result.at("options").at("dynamics").string() == "periodic-forcing";
    if (periodic) figure.references.erase(std::remove_if(figure.references.begin(), figure.references.end(),
        [](const auto &r) { return r.equilibrium; }), figure.references.end());
    const auto n = static_cast<size_t>(result.at("state_dimension").number());
    if (n < 2 || n > 64)
        return;
    const auto native =
        displayCoordinateTransform(metadata, n, result.at("coordinate_system").string());
    const auto display = displayCoordinateTransform(metadata, n, figure.coordinates);
    if (native.coordinates != result.at("coordinate_system").string())
        return;
    const double sx = display.stateMultipliers[0] / native.stateMultipliers[0],
                 sy = display.stateMultipliers[1] / native.stateMultipliers[1];
    auto add = [&](Trace trace) {
        if (trace.x.empty())
            return;
        if (figure.alternateCoordinates)
            figure.alternateCoordinates->traceLabels.push_back(trace.label);
        figure.traces.push_back(std::move(trace));
    };
    for (const auto &p : result.at("fixed_points").array()) {
        const auto &state = p.at("state").array();
        if (state.size() != n || !finite(state[0]) || !finite(state[1]))
            continue;
        const double x = state[0].number() * sx, y = state[1].number() * sy;
        const auto &lo = result.at("options").at("lower").array();
        const auto &hi = result.at("options").at("upper").array();
        const double tolerance = lo.size() == n && hi.size() == n
                                     ? 1e-7 * std::max((hi[0].number() - lo[0].number()) * sx,
                                                       (hi[1].number() - lo[1].number()) * sy)
                                     : 1e-9;
        if (std::none_of(figure.references.begin(), figure.references.end(), [&](const auto &r) {
                return r.equilibrium && std::hypot(r.x - x, r.y - y) < tolerance;
            }))
            figure.references.push_back({x, y, "EC", true});
    }
    auto curve = [&](const json::Value &c, const char *role, const char *label, const char *color) {
        if (c.at("status").string() == "unresolved")
            return;
        Trace trace;
        trace.role = role;
        trace.label = label;
        trace.color = color;
        trace.lineStyle = "dash";
        const auto &data = c.at("values").array();
        if (data.size() % n)
            return;
        for (size_t i = 0; i < data.size(); i += n)
            if (finite(data[i]) && finite(data[i + 1])) {
                trace.x.push_back(data[i].number() * sx);
                trace.y.push_back(data[i + 1].number() * sy);
            }
        for (auto &old : figure.traces)
            if (old.role == trace.role && old.color == trace.color && old.label == trace.label) {
                const double gap = std::numeric_limits<double>::quiet_NaN();
                old.x.push_back(gap);
                old.y.push_back(gap);
                old.x.insert(old.x.end(), trace.x.begin(), trace.x.end());
                old.y.insert(old.y.end(), trace.y.begin(), trace.y.end());
                return;
            }
        add(std::move(trace));
    };
    for (const auto &m : result.at("manifolds").array()) {
        const bool stable = m.at("kind").string() == "stable_manifold" || m.at("kind").string() == "stable_section_manifold";
        curve(m, stable ? "phase_stable_manifold" : "phase_unstable_manifold",
              periodic ? (stable ? "stable manifold (forcing section)" : "unstable manifold (forcing section)")
                       : (stable ? "stable manifold" : "unstable manifold"), stable ? "#668994" : "#ab8970");
    }
    appendPhaseSurfaces(figure, result, n, sx, sy);
    for (const auto &p : result.at("periodic_orbits").array())
        curve(p.at("curve"), "phase_periodic_orbit",
              periodic ? "periodic response" : (p.at("isolated").boolean() ? "limit cycle" : "periodic orbit"), "#859378");
    const auto &field = result.at("vector_field");
    const auto &axes = field.at("axes").array();
    const auto &points = field.at("points").array();
    const auto &vectors = field.at("vectors").array();
    const auto &low = result.at("options").at("lower").array();
    const auto &high = result.at("options").at("upper").array();
    if (axes.size() != 2 || axes[0].number() != 0 || axes[1].number() != 1 ||
        points.size() != vectors.size() || points.size() % n || low.size() != n || high.size() != n)
        return;
    const auto visible = automaticAxisLimits(figure);
    const double width = std::min((high[0].number() - low[0].number()) * sx, visible.xMax - visible.xMin),
                 height = std::min((high[1].number() - low[1].number()) * sy, visible.yMax - visible.yMin);
    if (!(width > 0 && height > 0)) return;
    const double length = .45 / result.at("options").at("vector_grid_per_axis").number(11);
    Trace arrows;
    arrows.role = "phase_vector_field";
    arrows.color = "#c1c8cc";
    for (size_t i = 0; i < points.size(); i += n) {
        if (!finite(points[i]) || !finite(points[i + 1]) || !finite(vectors[i]) ||
            !finite(vectors[i + 1]))
            continue;
        double dx = vectors[i].number() * sx / width, dy = vectors[i + 1].number() * sy / height;
        const double speed = std::hypot(dx, dy);
        if (speed < 1e-15)
            continue;
        dx *= length / speed;
        dy *= length / speed;
        const double x = points[i].number() * sx, y = points[i + 1].number() * sy,
                     u = x + dx * width, v = y + dy * height;
        segment(arrows, x, y, u, v);
        segment(arrows, u, v, u - (.3 * dx + .16 * dy) * width, v - (.3 * dy - .16 * dx) * height);
        segment(arrows, u, v, u - (.3 * dx - .16 * dy) * width, v - (.3 * dy + .16 * dx) * height);
    }
    add(std::move(arrows));
}
} // namespace dsl
