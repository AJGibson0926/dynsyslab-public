#include <dynsyslab/core/PathCrossings.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace dynsyslab::core {
namespace {
double det(double x, double y, double u, double v) { return x * v - y * u; }

struct Box { double x0, y0, x1, y1; };

struct Search {
    const std::vector<double>& x;
    const std::vector<double>& y;
    std::vector<PathCrossing>& out;
    std::size_t n;

    Box bounds(std::size_t i) const {
        return {std::min(x[i], x[i + 1]), std::min(y[i], y[i + 1]), std::max(x[i], x[i + 1]), std::max(y[i], y[i + 1])};
    }
    void test(std::size_t a, std::size_t b) const {
        const auto i = std::min(a, b), j = std::max(a, b);
        if (j - i <= 1) return;
        const double px = x[i + 1] - x[i], py = y[i + 1] - y[i];
        const double qx = x[j + 1] - x[j], qy = y[j + 1] - y[j];
        const double d = det(px, py, qx, qy);
        if (std::abs(d) <= 32 * std::numeric_limits<double>::epsilon() * std::hypot(px, py) * std::hypot(qx, qy)) return;
        const double dx = x[j] - x[i], dy = y[j] - y[i];
        const double f = det(dx, dy, qx, qy) / d, g = det(dx, dy, px, py) / d;
        if (f < 0 || g < 0 || f >= 1 || g > 1 || (g == 1 && j + 2 < n)) return;
        out.push_back({i, j, f, g});
    }
    // Adaptive quadtree: split crowded cells until each holds a few segments.
    // Dense regions (a decaying spiral near its sink) and long segments (a
    // collapse spike) both stay cheap: a segment only descends into the
    // quadrants its bounding box touches.
    void visit(const Box& cell, std::vector<std::size_t>& ids, int depth) const {
        if (ids.size() < 2) return;
        if (ids.size() <= 32 || depth >= 40) {
            for (std::size_t a = 0; a < ids.size(); ++a)
                for (std::size_t b = a + 1; b < ids.size(); ++b) {
                    const auto p = bounds(ids[a]), q = bounds(ids[b]);
                    if (p.x1 < q.x0 || q.x1 < p.x0 || p.y1 < q.y0 || q.y1 < p.y0) continue;
                    test(ids[a], ids[b]);
                }
            return;
        }
        const double mx = 0.5 * (cell.x0 + cell.x1), my = 0.5 * (cell.y0 + cell.y1);
        const Box quads[4]{{cell.x0, cell.y0, mx, my}, {mx, cell.y0, cell.x1, my},
                           {cell.x0, my, mx, cell.y1}, {mx, my, cell.x1, cell.y1}};
        std::vector<std::size_t> child[4];
        std::size_t total = 0;
        for (const auto i : ids) {
            const auto b = bounds(i);
            for (int k = 0; k < 4; ++k)
                if (b.x1 >= quads[k].x0 && b.x0 <= quads[k].x1 && b.y1 >= quads[k].y0 && b.y0 <= quads[k].y1) {
                    child[k].push_back(i); ++total; }
        }
        // Splitting no longer separates these segments (they all span the cell).
        if (total >= 4 * ids.size()) { visit(cell, ids, 40); return; }
        std::vector<std::size_t>().swap(ids);
        for (int k = 0; k < 4; ++k) visit(quads[k], child[k], depth + 1);
    }
};
}  // namespace

std::vector<PathCrossing> findPathCrossings(const std::vector<double>& x, const std::vector<double>& y) {
    const std::size_t n = x.size();
    if (n != y.size()) throw std::invalid_argument("crossing search needs matching x and y");
    std::vector<PathCrossing> out;
    if (n < 4) return out;
    for (std::size_t k = 0; k < n; ++k)
        if (!std::isfinite(x[k]) || !std::isfinite(y[k])) throw std::invalid_argument("crossing search needs finite points");
    std::vector<std::size_t> ids;
    ids.reserve(n - 1);
    for (std::size_t i = 0; i + 1 < n; ++i)
        if (x[i] != x[i + 1] || y[i] != y[i + 1]) ids.push_back(i);
    const Box root{*std::min_element(x.begin(), x.end()), *std::min_element(y.begin(), y.end()),
                   *std::max_element(x.begin(), x.end()), *std::max_element(y.begin(), y.end())};
    Search{x, y, out, n}.visit(root, ids, 0);
    // A pair that shares several leaves is found more than once.
    std::sort(out.begin(), out.end(), [](const PathCrossing& a, const PathCrossing& b) {
        return a.first != b.first ? a.first < b.first : a.second < b.second;
    });
    out.erase(std::unique(out.begin(), out.end(), [](const PathCrossing& a, const PathCrossing& b) {
        return a.first == b.first && a.second == b.second;
    }), out.end());
    return out;
}

} // namespace dynsyslab::core
