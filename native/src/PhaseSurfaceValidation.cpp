#include "PhaseSurfaceInternal.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace dynsyslab::core::phase_detail::surface {
namespace {
PhaseState vertex(const PhaseSurface &mesh, std::size_t i, std::size_t n) {
    return {mesh.vertices.begin() + i * n, mesh.vertices.begin() + (i + 1) * n};
}
double segmentDistance(const PhaseState &p, const PhaseState &a, const PhaseState &b,
                       const PhaseSpaceOptions &o) {
    auto edge = subtract(b, a), delta = subtract(p, a);
    const double square = dot(edge, edge, o);
    const double t = square > 0 ? std::clamp(dot(delta, edge, o) / square, 0., 1.) : 0.;
    for (std::size_t j = 0; j < p.size(); ++j)
        delta[j] -= t * edge[j];
    return length(delta, o);
}
double triangleDistance(const PhaseState &p, const PhaseState &a, const PhaseState &b,
                        const PhaseState &c, const PhaseSpaceOptions &o) {
    const auto u = subtract(b, a), v = subtract(c, a), d = subtract(p, a);
    const double aa = dot(u, u, o), bb = dot(v, v, o), ab = dot(u, v, o), det = aa * bb - ab * ab;
    if (det > 1e-12 * aa * bb) {
        const double du = dot(d, u, o), dv = dot(d, v, o), s = (bb * du - ab * dv) / det,
                     t = (aa * dv - ab * du) / det;
        if (s >= 0 && t >= 0 && s + t <= 1) {
            auto residual = d;
            for (std::size_t j = 0; j < p.size(); ++j)
                residual[j] -= s * u[j] + t * v[j];
            return length(residual, o);
        }
    }
    return std::min(
        {segmentDistance(p, a, b, o), segmentDistance(p, b, c, o), segmentDistance(p, c, a, o)});
}
double nearest(const PhaseState &point, const PhaseSurface &mesh, const PhaseSpaceOptions &o) {
    const auto n = point.size();
    double distance = std::numeric_limits<double>::infinity();
    for (std::size_t k = 0; k < mesh.triangles.size(); k += 3)
        distance = std::min(distance, triangleDistance(point, vertex(mesh, mesh.triangles[k], n),
                                                       vertex(mesh, mesh.triangles[k + 1], n),
                                                       vertex(mesh, mesh.triangles[k + 2], n), o));
    return distance;
}
} // namespace
void validate(const ODESystem &system, PhaseSurface &mesh, const PhaseSurface &refined,
              const PhaseSpaceOptions &o, const ExecutionContext *execution) {
    const auto n = o.lower.size(), count = mesh.vertices.size() / n,
               faces = mesh.triangles.size() / 3;
    if (mesh.frontOffsets.size() < 3 || refined.frontOffsets.size() < 3) {
        mesh.status = "unresolved";
        mesh.reason += "; insufficient continued fronts for independent validation";
        return;
    }
    // Full-state field tangency is checked at every face centroid, independently of shooting.
    for (std::size_t k = 0; k < faces; ++k) {
        if (execution)
            execution->checkpoint();
        const auto a = vertex(mesh, mesh.triangles[3 * k], n),
                   b = vertex(mesh, mesh.triangles[3 * k + 1], n),
                   c = vertex(mesh, mesh.triangles[3 * k + 2], n);
        auto center = a;
        for (std::size_t j = 0; j < n; ++j)
            center[j] = (a[j] + b[j] + c[j]) / 3;
        const auto u = subtract(b, a), v = subtract(c, a);
        auto f = field(system, center);
        const double aa = dot(u, u, o), bb = dot(v, v, o), ab = dot(u, v, o),
                     det = aa * bb - ab * ab;
        if (!std::isfinite(det) || det <= 1e-12 * aa * bb) {
            mesh.status = "unresolved";
            mesh.reason += "; mesh tangent chart is rank-deficient";
            return;
        }
        const double fu = dot(f, u, o), fv = dot(f, v, o), s = (bb * fu - ab * fv) / det,
                     t = (aa * fv - ab * fu) / det;
        const double speed = length(f, o);
        for (std::size_t j = 0; j < n; ++j)
            f[j] -= s * u[j] + t * v[j];
        mesh.tangentDefect = std::max(mesh.tangentDefect, length(f, o) / std::max(speed, 1e-12));
        if (k % std::max(std::size_t(1), faces / 256) == 0)
            mesh.refinementError = std::max(mesh.refinementError, nearest(center, refined, o));
    }
    // Include the complete boundary and distributed interior vertices in the geometric check.
    const auto boundary = mesh.frontOffsets.back();
    for (std::size_t k = 0; k < count; ++k) {
        if (k < boundary && k % std::max(std::size_t(1), count / 256) != 0)
            continue;
        if (execution)
            execution->checkpoint();
        mesh.refinementError =
            std::max(mesh.refinementError, nearest(vertex(mesh, k, n), refined, o));
    }
    if (!std::isfinite(mesh.refinementError) ||
        mesh.refinementError > 4 * o.surface.chordTolerance ||
        mesh.tangentDefect > o.surface.tangentTolerance) {
        mesh.status = "unresolved";
        mesh.reason += "; seed/mesh refinement or off-sample field tangency did not converge";
    }
}
} // namespace dynsyslab::core::phase_detail::surface
