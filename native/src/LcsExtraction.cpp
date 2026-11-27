// FTLE ridges and variational LCS (strainlines) on two-axis grids (LagrangianCoherence.hpp).
//
// Ridges: Shadden, Lekien and Marsden 2005, Physica D 212:271-304, section 4 (second-derivative
// ridge: a maximum of the field along the Hessian eigenvector of most negative curvature).
// Strainlines: Haller 2011, Physica D 240:574-598 (variational LCS conditions) and Farazmand and
// Haller 2012, Chaos 22:013128 (strainlines r' = xi_1, conditions A, B and D, failure length l_f,
// and selection of the segments with locally maximal mean lambda_2).
#include <dynsyslab/core/LagrangianCoherence.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

struct Plane {
    std::size_t nx{}, ny{};
    double x0{}, y0{}, hx{}, hy{};
    explicit Plane(const DomainGrid& grid) {
        if (grid.axes.size() != 2) throw std::invalid_argument("LCS extraction needs a grid of two axes");
        nx = grid.axes[0].points; ny = grid.axes[1].points;
        if (nx < 3 || ny < 3) throw std::invalid_argument("LCS extraction needs at least three nodes per axis");
        x0 = grid.axes[0].lower; y0 = grid.axes[1].lower;
        hx = (grid.axes[0].upper - x0) / double(nx - 1);
        hy = (grid.axes[1].upper - y0) / double(ny - 1);
    }
    [[nodiscard]] bool inside(double x, double y) const {
        return x >= x0 && y >= y0 && x <= x0 + hx * double(nx - 1) && y <= y0 + hy * double(ny - 1);
    }
    [[nodiscard]] std::size_t nearest(double x, double y) const {
        const auto i = std::size_t(std::clamp(std::lround((x - x0) / hx), 0L, long(nx - 1)));
        const auto j = std::size_t(std::clamp(std::lround((y - y0) / hy), 0L, long(ny - 1)));
        return i + nx * j;
    }
    // Bilinear interpolation of node values f (NaN when a corner is NaN).
    [[nodiscard]] double sample(const std::vector<double>& f, double x, double y) const {
        const double u = std::clamp((x - x0) / hx, 0.0, double(nx - 1)), v = std::clamp((y - y0) / hy, 0.0, double(ny - 1));
        const auto i = std::min(std::size_t(u), nx - 2), j = std::min(std::size_t(v), ny - 2);
        const double a = u - double(i), b = v - double(j);
        return (1 - a) * (1 - b) * f[i + nx * j] + a * (1 - b) * f[i + 1 + nx * j]
             + (1 - a) * b * f[i + nx * (j + 1)] + a * b * f[i + 1 + nx * (j + 1)];
    }
    // Gradient and Hessian of f by central differences at an interior node; false if any value is NaN.
    bool derivatives(const std::vector<double>& f, std::size_t i, std::size_t j, double& fx, double& fy,
                     double& fxx, double& fyy, double& fxy) const {
        const auto at = [&](std::size_t a, std::size_t b) { return f[a + nx * b]; };
        const double c = at(i, j), e = at(i + 1, j), w = at(i - 1, j), n = at(i, j + 1), s = at(i, j - 1);
        const double ne = at(i + 1, j + 1), nw = at(i - 1, j + 1), se = at(i + 1, j - 1), sw = at(i - 1, j - 1);
        fx = (e - w) / (2 * hx); fy = (n - s) / (2 * hy);
        fxx = (e - 2 * c + w) / (hx * hx); fyy = (n - 2 * c + s) / (hy * hy);
        fxy = (ne - nw - se + sw) / (4 * hx * hy);
        return std::isfinite(fx + fy + fxx + fyy + fxy);
    }
};

// Most negative curvature direction of the 2 x 2 Hessian; false when its curvature is not negative.
bool ridgeNormal(double fxx, double fyy, double fxy, double& nx, double& ny) {
    std::vector<double> values, vectors;
    symmetricEigen({fxx, fxy, fxy, fyy}, 2, values, vectors);
    nx = vectors[0]; ny = vectors[2];
    return values[0] < 0;
}

}  // namespace

std::vector<std::uint8_t> ridgeNodes(const DomainGrid& grid, const std::vector<double>& f, double threshold) {
    const Plane plane(grid);
    if (f.size() != plane.nx * plane.ny) throw std::invalid_argument("ridgeNodes: one value per grid node");
    double top = -std::numeric_limits<double>::infinity();
    for (const auto v : f) if (std::isfinite(v)) top = std::max(top, v);
    std::vector<std::uint8_t> out(f.size(), 0);
    if (!std::isfinite(top)) return out;
    const double floor = threshold * top, h = std::min(plane.hx, plane.hy);
    for (std::size_t j = 1; j + 1 < plane.ny; ++j)
        for (std::size_t i = 1; i + 1 < plane.nx; ++i) {
            const double c = f[i + plane.nx * j];
            double fx, fy, fxx, fyy, fxy, nx, ny;
            if (!std::isfinite(c) || c < floor || !plane.derivatives(f, i, j, fx, fy, fxx, fyy, fxy)) continue;
            if (!ridgeNormal(fxx, fyy, fxy, nx, ny)) continue;
            const double x = plane.x0 + plane.hx * double(i), y = plane.y0 + plane.hy * double(j);
            // A maximum along the normal: ties broken to one side so a ridge between nodes marks one.
            const double ahead = plane.sample(f, x + h * nx, y + h * ny), behind = plane.sample(f, x - h * nx, y - h * ny);
            if (c > ahead && c >= behind) out[i + plane.nx * j] = 1;
        }
    return out;
}

std::vector<LagrangianCurve> strainlines(const CauchyGreenField& field, const StrainlineOptions& options) {
    const Plane plane(field.grid);
    if (field.axes != 2) throw std::invalid_argument("strainlines need a two-axis Cauchy-Green field");
    const auto nodes = plane.nx * plane.ny;
    std::vector<double> l1(nodes), l2(nodes), xi1x(nodes), xi1y(nodes);
    for (std::size_t p = 0; p < nodes; ++p) {
        l1[p] = field.eigenvalues[p * 2];
        l2[p] = field.eigenvalues[p * 2 + 1];
        xi1x[p] = field.eigenvectors[p * 4 + 0];  // column 0: (row 0, row 1)
        xi1y[p] = field.eigenvectors[p * 4 + 2];
    }
    // Conditions A (lambda_1 != lambda_2 > 1) and D (<xi_2, D^2 lambda_2 xi_2> < 0) at the nodes.
    std::vector<std::uint8_t> admissible(nodes, 0);
    for (std::size_t j = 1; j + 1 < plane.ny; ++j)
        for (std::size_t i = 1; i + 1 < plane.nx; ++i) {
            const auto p = i + plane.nx * j;
            double fx, fy, fxx, fyy, fxy;
            if (!(l2[p] > 1) || !(l2[p] - l1[p] > 1e-9 * l2[p]) || !plane.derivatives(l2, i, j, fx, fy, fxx, fyy, fxy)) continue;
            const double ax = -xi1y[p], ay = xi1x[p];  // xi_2 is orthogonal to xi_1
            admissible[p] = (ax * ax * fxx + 2 * ax * ay * fxy + ay * ay * fyy) < 0;
        }
    const double h = std::min(plane.hx, plane.hy);
    const double width = plane.hx * double(plane.nx - 1), height = plane.hy * double(plane.ny - 1);
    const double ds = options.step * h;
    const double maximum = options.maximumLength > 0 ? options.maximumLength : std::hypot(width, height);
    const double minimum = options.minimumLength > 0 ? options.minimumLength : std::min(width, height) / 5;
    const double failure = options.failureLength > 0 ? options.failureLength : 2 * h;
    const double separation = options.separation > 0 ? options.separation : 2 * h;
    if (!(ds > 0)) throw std::invalid_argument("strainlines: step must be positive");

    // xi_1 at (x, y), its sign aligned with the direction travelled (corners aligned first).
    const auto direction = [&](double x, double y, double dx, double dy, double& ox, double& oy) {
        const double u = std::clamp((x - plane.x0) / plane.hx, 0.0, double(plane.nx - 1));
        const double v = std::clamp((y - plane.y0) / plane.hy, 0.0, double(plane.ny - 1));
        const auto i = std::min(std::size_t(u), plane.nx - 2), j = std::min(std::size_t(v), plane.ny - 2);
        const double a = u - double(i), b = v - double(j);
        ox = oy = 0;
        const std::size_t corner[4] = {i + plane.nx * j, i + 1 + plane.nx * j, i + plane.nx * (j + 1), i + 1 + plane.nx * (j + 1)};
        const double weight[4] = {(1 - a) * (1 - b), a * (1 - b), (1 - a) * b, a * b};
        for (int c = 0; c < 4; ++c) {
            double ex = xi1x[corner[c]], ey = xi1y[corner[c]];
            if (!std::isfinite(ex + ey)) return false;
            if (ex * dx + ey * dy < 0) { ex = -ex; ey = -ey; }
            ox += weight[c] * ex; oy += weight[c] * ey;
        }
        const double norm = std::hypot(ox, oy);
        if (!(norm > 0)) return false;
        ox /= norm; oy /= norm;
        return true;
    };
    // One branch from the seed along +-xi_1: samples and whether each satisfies the conditions.
    const auto branch = [&](double x, double y, double dx, double dy, std::vector<double>& pts, std::vector<bool>& ok) {
        double travelled = 0, failing = 0;
        while (travelled < maximum / 2) {
            double k1x, k1y, k2x, k2y, k3x, k3y, k4x, k4y;
            if (!direction(x, y, dx, dy, k1x, k1y) || !direction(x + ds / 2 * k1x, y + ds / 2 * k1y, k1x, k1y, k2x, k2y)
                || !direction(x + ds / 2 * k2x, y + ds / 2 * k2y, k1x, k1y, k3x, k3y)
                || !direction(x + ds * k3x, y + ds * k3y, k1x, k1y, k4x, k4y)) break;
            const double sx = (k1x + 2 * k2x + 2 * k3x + k4x) / 6, sy = (k1y + 2 * k2y + 2 * k3y + k4y) / 6;
            x += ds * sx; y += ds * sy;
            if (!plane.inside(x, y)) break;
            dx = sx; dy = sy;
            travelled += ds;
            const bool good = admissible[plane.nearest(x, y)] != 0;
            failing = good ? 0 : failing + ds;
            if (failing > failure) break;
            pts.push_back(x); pts.push_back(y); ok.push_back(good);
        }
    };

    const auto seeds = ridgeNodes(field.grid, l2, options.threshold);
    std::vector<std::size_t> order;
    for (std::size_t p = 0; p < nodes; ++p) {
        const auto i = p % plane.nx, j = p / plane.nx;
        if (seeds[p] && admissible[p] && (i + j) % std::max<std::size_t>(1, options.seedStride) == 0) order.push_back(p);
    }
    std::sort(order.begin(), order.end(), [&](auto a, auto b) { return l2[a] > l2[b]; });

    std::vector<LagrangianCurve> candidates;
    std::vector<std::uint8_t> visited(nodes, 0);
    for (const auto p : order) {
        if (visited[p]) continue;  // already on a computed strainline
        const double x = plane.x0 + plane.hx * double(p % plane.nx), y = plane.y0 + plane.hy * double(p / plane.nx);
        std::vector<double> forward, backward;
        std::vector<bool> okForward, okBackward;
        branch(x, y, xi1x[p], xi1y[p], forward, okForward);
        branch(x, y, -xi1x[p], -xi1y[p], backward, okBackward);
        // Join: reversed backward branch, the seed, forward branch; trim failing ends.
        std::vector<double> pts;
        std::vector<bool> ok;
        for (std::size_t s = okBackward.size(); s-- > 0;) {
            pts.push_back(backward[2 * s]); pts.push_back(backward[2 * s + 1]); ok.push_back(okBackward[s]);
        }
        pts.push_back(x); pts.push_back(y); ok.push_back(true);
        pts.insert(pts.end(), forward.begin(), forward.end());
        ok.insert(ok.end(), okForward.begin(), okForward.end());
        std::size_t first = 0, last = ok.size();
        while (first < last && !ok[first]) ++first;
        while (last > first && !ok[last - 1]) --last;
        LagrangianCurve curve;
        curve.points.assign(pts.begin() + std::ptrdiff_t(2 * first), pts.begin() + std::ptrdiff_t(2 * last));
        const auto samples = curve.points.size() / 2;
        double total = 0;
        for (std::size_t s = 0; s < samples; ++s) {
            const double cx = curve.points[2 * s], cy = curve.points[2 * s + 1];
            total += plane.sample(l2, cx, cy);
            visited[plane.nearest(cx, cy)] = 1;
            if (s) curve.length += std::hypot(cx - curve.points[2 * s - 2], cy - curve.points[2 * s - 1]);
        }
        if (samples < 2 || curve.length < minimum) continue;
        curve.meanStretching = total / double(samples);
        candidates.push_back(std::move(curve));
    }
    // Keep the strongest of nearby candidates: a curve mostly within `separation` of a kept one goes.
    std::sort(candidates.begin(), candidates.end(), [](const auto& a, const auto& b) { return a.meanStretching > b.meanStretching; });
    const auto rx = long(std::ceil(separation / plane.hx)), ry = long(std::ceil(separation / plane.hy));
    std::vector<std::uint8_t> taken(nodes, 0);
    std::vector<LagrangianCurve> out;
    for (auto& curve : candidates) {
        const auto samples = curve.points.size() / 2;
        std::size_t near = 0;
        for (std::size_t s = 0; s < samples; ++s) near += taken[plane.nearest(curve.points[2 * s], curve.points[2 * s + 1])];
        if (2 * near > samples) continue;
        for (std::size_t s = 0; s < samples; ++s) {
            const auto q = plane.nearest(curve.points[2 * s], curve.points[2 * s + 1]);
            const long i = long(q % plane.nx), j = long(q / plane.nx);
            for (long b = std::max(0L, j - ry); b <= std::min(long(plane.ny) - 1, j + ry); ++b)
                for (long a = std::max(0L, i - rx); a <= std::min(long(plane.nx) - 1, i + rx); ++a)
                    taken[std::size_t(a) + plane.nx * std::size_t(b)] = 1;
        }
        out.push_back(std::move(curve));
    }
    return out;
}

std::vector<std::uint8_t> rasterizeCurves(const DomainGrid& grid, const std::vector<LagrangianCurve>& curves) {
    const Plane plane(grid);
    std::vector<std::uint8_t> out(plane.nx * plane.ny, 0);
    for (const auto& curve : curves)
        for (std::size_t s = 0; s + 1 < curve.points.size(); s += 2) out[plane.nearest(curve.points[s], curve.points[s + 1])] = 1;
    return out;
}

}  // namespace dynsyslab::core
