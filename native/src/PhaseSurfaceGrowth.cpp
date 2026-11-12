#include "PhaseSurfaceInternal.hpp"
#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>
namespace dynsyslab::core::phase_detail::surface {
namespace {
void append(PhaseSurface &mesh, const PhaseState &x) {
    mesh.vertices.insert(mesh.vertices.end(), x.begin(), x.end());
}
void triangle(PhaseSurface &mesh, std::size_t a, std::size_t b, std::size_t c) {
    mesh.triangles.insert(mesh.triangles.end(), {a, b, c});
}
std::vector<Point> nextFront(const ODESystem &system, const Front &front, double step, double sign,
                             const PhaseSpaceOptions &o, Budget &budget,
                             const ExecutionContext *execution) {
    std::vector<Point> anchors, points;
    for (std::size_t j = 0; j < front.points.size(); ++j)
        anchors.push_back(advance(system, front, double(j), step, sign, o, budget, execution));
    std::function<void(const Point &, const Point &, int)> edge;
    edge = [&](const Point &a, const Point &b, int depth) {
        if (points.size() >= o.surface.maximumPoints)
            throw std::runtime_error("Surface front point budget reached");
        auto midpoint =
            advance(system, front, (a.parent + b.parent) / 2, step, sign, o, budget, execution);
        auto chord = a.state;
        for (std::size_t k = 0; k < chord.size(); ++k)
            chord[k] = (chord[k] + b.state[k]) / 2;
        if (length(subtract(a.state, b.state), o) > o.surface.meshSize ||
            length(subtract(midpoint.state, chord), o) > o.surface.chordTolerance) {
            if (depth >= 10)
                throw std::runtime_error("Surface front interpolation could not be resolved");
            edge(a, midpoint, depth + 1);
            edge(midpoint, b, depth + 1);
        } else
            points.push_back(a);
    };
    for (std::size_t j = 0; j < anchors.size(); ++j) {
        auto b = anchors[(j + 1) % anchors.size()];
        b.parent = double(j + 1);
        edge(anchors[j], b, 0);
    }
    return points;
}
void ribbon(PhaseSurface &mesh, const std::vector<Point> &next, std::size_t oldStart,
            std::size_t oldCount, std::size_t n) {
    const auto start = mesh.vertices.size() / n;
    for (const auto &p : next)
        append(mesh, p.state);
    std::size_t cursor = 0;
    for (std::size_t j = 0; j < oldCount; ++j) {
        // Insertion-only refinement keeps each old integer parameter in the new ring.
        while (cursor + 1 < next.size() && next[cursor + 1].parent <= double(j + 1)) {
            triangle(mesh, oldStart + j, start + cursor, start + cursor + 1);
            ++cursor;
        }
        if (j + 1 == oldCount) {
            triangle(mesh, oldStart + j, start + cursor, start);
            triangle(mesh, oldStart + j, start, oldStart);
        } else
            triangle(mesh, oldStart + j, start + cursor, oldStart + j + 1);
    }
    mesh.frontOffsets.push_back(start);
}
} // namespace
PhaseSurface grow(const ODESystem &system, const PhaseState &center, const Ring &basis, double sign,
                  const PhaseSpaceOptions &o, const ExecutionContext *execution) {
    PhaseSurface mesh;
    mesh.status = "truncated";
    mesh.reason = "Surface front budget reached";
    mesh.seedRadius = o.surface.seedRadius;
    const auto n = center.size();
    Front front;
    Budget budget;
    append(mesh, center);
    mesh.frontOffsets.push_back(1);
    mesh.frontDistances.push_back(o.surface.seedRadius);
    for (std::size_t j = 0; j < o.surface.initialPoints; ++j) {
        const double theta = 2 * std::acos(-1.) * j / o.surface.initialPoints;
        auto x = center;
        for (std::size_t k = 0; k < n; ++k)
            x[k] += o.surface.seedRadius *
                    (std::cos(theta) * basis[0][k] + std::sin(theta) * basis[1][k]);
        if (!inside(x, o))
            throw std::runtime_error("Surface seed ellipse leaves the declared domain");
        front.points.push_back(x);
        front.previous.push_back(center);
        append(mesh, x);
        triangle(mesh, 0, 1 + j, 1 + (j + 1) % o.surface.initialPoints);
    }
    double distance = o.surface.seedRadius,
           step = std::min(o.surface.step, o.surface.seedRadius / 2);
    for (std::size_t ring = 1; ring < o.surface.maximumFronts; ++ring) {
        if (execution)
            execution->checkpoint();
        step = std::min(step, o.surface.maximumDistance - distance);
        if (step <= 1e-10) {
            mesh.status = "resolved";
            mesh.reason = "Requested finite surface extent reached";
            break;
        }
        bool accepted = false;
        std::vector<Point> next;
        double maximumTurn = 0;
        for (int trial = 0; trial < 9 && !accepted; ++trial) {
            try {
                next = nextFront(system, front, step, sign, o, budget, execution);
                maximumTurn = 0;
                for (const auto &p : next) {
                    const auto r = interpolate(front.points, p.parent);
                    maximumTurn = std::max(maximumTurn,
                                           angle(subtract(r, interpolate(front.previous, p.parent)),
                                                 subtract(p.state, r), o));
                }
                // The seed ellipse has no preceding geodesic increment; its radius is not a normal.
                if (ring == 1)
                    maximumTurn = 0;
                if (maximumTurn > .55 || step * maximumTurn > 4 * o.surface.chordTolerance)
                    throw std::runtime_error(
                        "Surface radial curvature exceeds the continuation tolerance");
                if (mesh.vertices.size() / n + next.size() > o.surface.maximumPoints) {
                    mesh.reason = "Surface vertex budget reached";
                    return mesh;
                }
                accepted = true;
            } catch (const ExecutionCancelled &) {
                throw;
            } catch (const std::exception &error) {
                if (trial == 0)
                    mesh.reason = error.what();
                step *= .5;
                if (step < o.surface.step * 1e-4)
                    break;
            }
        }
        if (!accepted)
            return mesh;
        Front updated;
        for (const auto &p : next) {
            updated.points.push_back(p.state);
            updated.previous.push_back(interpolate(front.points, p.parent));
            mesh.shootingResidual = std::max(mesh.shootingResidual, p.residual);
        }
        ribbon(mesh, next, mesh.frontOffsets.back(), front.points.size(), n);
        distance += step;
        mesh.frontDistances.push_back(distance);
        mesh.maximumTurn = std::max(mesh.maximumTurn, maximumTurn);
        front = std::move(updated);
        mesh.reason = "Surface front budget reached";
        if (maximumTurn < .15 && step * maximumTurn < o.surface.chordTolerance)
            step = std::min(o.surface.step, step * 1.5);
    }
    if (distance >= o.surface.maximumDistance - 1e-10) {
        mesh.status = "resolved";
        mesh.reason = "Requested finite surface extent reached";
    }
    return mesh;
}
} // namespace dynsyslab::core::phase_detail::surface
