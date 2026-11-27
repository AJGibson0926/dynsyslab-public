// Planar orbit geometry for injectivity audits (OrbitGeometry.hpp).
#include <dynsyslab/core/OrbitGeometry.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace dynsyslab::core {

std::vector<std::size_t> findPeaks(const std::vector<double>& y, std::size_t distance, double prominence) {
    const auto n = y.size();
    std::vector<std::size_t> candidates;
    for (std::size_t i = 1; i + 1 < n; ++i) {
        if (!(y[i] > y[i - 1])) continue;
        std::size_t j = i;
        while (j + 1 < n && y[j + 1] == y[i]) ++j;  // a plateau counts once
        if (j + 1 < n && y[j + 1] < y[i]) candidates.push_back(i);
        i = j;
    }
    std::vector<std::size_t> prominent;
    for (const auto p : candidates) {
        double left = y[p], right = y[p];
        for (std::size_t k = p; k-- > 0 && y[k] <= y[p];) left = std::min(left, y[k]);
        for (std::size_t k = p + 1; k < n && y[k] <= y[p]; ++k) right = std::min(right, y[k]);
        if (y[p] - std::max(left, right) >= prominence) prominent.push_back(p);
    }
    if (distance <= 1) return prominent;
    std::vector<std::size_t> order(prominent.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](auto a, auto b) { return y[prominent[a]] > y[prominent[b]]; });
    std::vector<bool> kept(prominent.size(), true);
    for (const auto o : order) {
        if (!kept[o]) continue;
        for (std::size_t k = 0; k < prominent.size(); ++k)
            if (k != o && kept[k] && (prominent[k] > prominent[o] ? prominent[k] - prominent[o] : prominent[o] - prominent[k]) < distance)
                kept[k] = false;
    }
    std::vector<std::size_t> out;
    for (std::size_t k = 0; k < prominent.size(); ++k) if (kept[k]) out.push_back(prominent[k]);
    return out;
}

std::vector<double> representativeCycle(const std::vector<double>& path, std::size_t width, const std::vector<std::size_t>& peaks,
                                        std::size_t samples, std::size_t which) {
    const auto rows = width ? path.size() / width : 0;
    if (peaks.size() < 3 + which || samples < 2) throw std::invalid_argument("representative cycle: too few complete cycles");
    const auto start = peaks[peaks.size() - 3 - which], stop = peaks[peaks.size() - 2 - which];
    if (stop <= start || stop >= rows) throw std::invalid_argument("representative cycle: peaks outside the path");
    std::vector<double> out(samples * width);
    for (std::size_t s = 0; s < samples; ++s) {
        const double at = double(start) + double(stop - start) * double(s) / double(samples - 1);
        const auto i = std::min(std::size_t(at), stop - 1);
        const double f = at - double(i);
        for (std::size_t c = 0; c < width; ++c) out[s * width + c] = path[i * width + c] + f * (path[(i + 1) * width + c] - path[i * width + c]);
    }
    return out;
}

double peakClosure(const std::vector<double>& path, std::size_t width, const std::vector<std::size_t>& peaks) {
    const auto rows = width ? path.size() / width : 0;
    if (peaks.size() < 3 || !rows) return std::numeric_limits<double>::quiet_NaN();
    std::vector<double> mean(width, 0);
    for (std::size_t k = 0; k < rows; ++k) for (std::size_t c = 0; c < width; ++c) mean[c] += path[k * width + c] / double(rows);
    double radius = 0;
    for (std::size_t k = 0; k < rows; ++k) for (std::size_t c = 0; c < width; ++c) radius += std::pow(path[k * width + c] - mean[c], 2);
    radius = std::sqrt(radius / double(rows));
    std::vector<double> gaps;
    for (std::size_t k = 1; k < peaks.size(); ++k) {
        double g = 0;
        for (std::size_t c = 0; c < width; ++c) g += std::pow(path[peaks[k] * width + c] - path[peaks[k - 1] * width + c], 2);
        gaps.push_back(std::sqrt(g));
    }
    std::sort(gaps.begin(), gaps.end());
    const auto m = gaps.size();
    const double median = m % 2 ? gaps[m / 2] : (gaps[m / 2 - 1] + gaps[m / 2]) / 2;
    return radius > 0 ? median / radius : std::numeric_limits<double>::quiet_NaN();
}

std::vector<CurveIntersection> findCurveIntersections(const std::vector<double>& ax, const std::vector<double>& ay,
                                                      const std::vector<double>& bx, const std::vector<double>& by, double merge) {
    if (ax.size() != ay.size() || bx.size() != by.size()) throw std::invalid_argument("curve intersections: x and y differ in length");
    std::vector<CurveIntersection> out;
    constexpr double edge = 1e-9;
    for (std::size_t i = 0; i + 1 < ax.size(); ++i) {
        const double px = ax[i], py = ay[i], rx = ax[i + 1] - px, ry = ay[i + 1] - py;
        const double xlo = std::min(px, px + rx), xhi = std::max(px, px + rx), ylo = std::min(py, py + ry), yhi = std::max(py, py + ry);
        for (std::size_t j = 0; j + 1 < bx.size(); ++j) {
            const double qx = bx[j], qy = by[j], sx = bx[j + 1] - qx, sy = by[j + 1] - qy;
            if (std::max(qx, qx + sx) < xlo || std::min(qx, qx + sx) > xhi || std::max(qy, qy + sy) < ylo || std::min(qy, qy + sy) > yhi) continue;
            const double d = rx * sy - ry * sx;
            if (std::abs(d) <= 1e-13) continue;  // parallel
            const double ox = qx - px, oy = qy - py;
            const double u = (ox * sy - oy * sx) / d, v = (ox * ry - oy * rx) / d;
            if (u <= edge || u >= 1 - edge || v <= edge || v >= 1 - edge) continue;
            const double x = px + u * rx, y = py + u * ry;
            const bool duplicate = std::any_of(out.begin(), out.end(), [&](const auto& c) { return std::hypot(c.x - x, c.y - y) < merge; });
            if (!duplicate) out.push_back({i, j, u, v, x, y});
        }
    }
    return out;
}

double polygonArea(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size()) throw std::invalid_argument("polygon area: x and y differ in length");
    double s = 0;
    for (std::size_t k = 0; k < x.size(); ++k) {
        const auto next = (k + 1) % x.size();
        s += x[k] * y[next] - y[k] * x[next];
    }
    return std::abs(s) / 2;
}

double nearestSeparation(const std::vector<double>& a, const std::vector<double>& b, std::size_t width) {
    if (!width || a.size() % width || b.size() % width) throw std::invalid_argument("nearest separation: path widths");
    double best = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < a.size(); i += width)
        for (std::size_t j = 0; j < b.size(); j += width) {
            double s = 0;
            for (std::size_t c = 0; c < width && s < best * best; ++c) s += std::pow(a[i + c] - b[j + c], 2);
            best = std::min(best, std::sqrt(s));
        }
    return best;
}

}  // namespace dynsyslab::core

namespace dynsyslab::core {
namespace {

// Value and slope of the linear interpolant of (time, v) at t (clamped to the record).
std::pair<double, double> interpolate(const std::vector<double>& time, const std::vector<double>& v, double t) {
    const auto j = std::min(std::max<std::size_t>(1, std::size_t(std::upper_bound(time.begin(), time.end(), t) - time.begin())), time.size() - 1);
    const double slope = (v[j] - v[j - 1]) / (time[j] - time[j - 1]);
    return {v[j - 1] + slope * (t - time[j - 1]), slope};
}

}  // namespace

RefinedVisit refineVisits(const std::vector<double>& time, const std::vector<double>& x, const std::vector<double>& y, double a0,
                          double a1, double b0, double b1, double sx, double sy, double tolerance) {
    if (time.size() < 2 || x.size() != time.size() || y.size() != time.size() || !(a0 < a1) || !(b0 < b1) || a1 > b0 || !(sx > 0) || !(sy > 0))
        throw std::invalid_argument("refine visits: a record and two ordered, distinct brackets");
    double t1 = (a0 + a1) / 2, t2 = (b0 + b1) / 2;
    RefinedVisit out;
    for (int iteration = 0; iteration < 60; ++iteration) {
        const auto [x1, dx1] = interpolate(time, x, t1);
        const auto [y1, dy1] = interpolate(time, y, t1);
        const auto [x2, dx2] = interpolate(time, x, t2);
        const auto [y2, dy2] = interpolate(time, y, t2);
        const double rx = (x1 - x2) / sx, ry = (y1 - y2) / sy;
        out.gap = std::hypot(rx, ry);
        if (out.gap <= tolerance) break;
        // J = [[dx1, -dx2], [dy1, -dy2]] / scales; Newton step J d = -r.
        const double j11 = dx1 / sx, j12 = -dx2 / sx, j21 = dy1 / sy, j22 = -dy2 / sy;
        const double det = j11 * j22 - j12 * j21;
        if (std::abs(det) < 1e-300) break;
        t1 = std::clamp(t1 + (-rx * j22 + ry * j12) / det, a0, a1);
        t2 = std::clamp(t2 + (-ry * j11 + rx * j21) / det, b0, b1);
    }
    out.first = t1; out.second = t2;
    out.matched = out.gap <= tolerance && t1 < t2;
    return out;
}

VisitMatching matchVisits(const std::vector<std::pair<double, double>>& a, const std::vector<std::pair<double, double>>& b,
                          double tol1, double tol2) {
    if (!(tol1 > 0) || !(tol2 > 0)) throw std::invalid_argument("match visits: positive tolerances");
    const auto m = a.size(), n = b.size();
    // Square cost matrix (m + n): admissible pairs cost their scaled shift, dummies a penalty larger
    // than any sum of admissible costs, inadmissible pairs twice that (as in the side-car).
    const double penalty = double(m + n + 1);
    const auto N = m + n;
    std::vector<double> cost(N * N, 0);
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = 0; j < N; ++j) {
            if (i < m && j < n) {
                const double c = std::max(std::abs(a[i].first - b[j].first) / tol1, std::abs(a[i].second - b[j].second) / tol2);
                cost[i * N + j] = c <= 1 ? c : 2 * penalty;
            } else if (i < m || j < n) {
                cost[i * N + j] = penalty;
            }
        }
    // Hungarian method with potentials (Kuhn 1955; O(N^3) formulation of Jonker and Volgenant 1987).
    const double inf = std::numeric_limits<double>::infinity();
    std::vector<double> u(N + 1, 0), v(N + 1, 0);
    std::vector<std::size_t> p(N + 1, 0), way(N + 1, 0);
    for (std::size_t i = 1; i <= N; ++i) {
        p[0] = i;
        std::size_t j0 = 0;
        std::vector<double> minv(N + 1, inf);
        std::vector<bool> used(N + 1, false);
        do {
            used[j0] = true;
            const auto i0 = p[j0];
            double delta = inf;
            std::size_t j1 = 0;
            for (std::size_t j = 1; j <= N; ++j) {
                if (used[j]) continue;
                const double cur = cost[(i0 - 1) * N + (j - 1)] - u[i0] - v[j];
                if (cur < minv[j]) { minv[j] = cur; way[j] = j0; }
                if (minv[j] < delta) { delta = minv[j]; j1 = j; }
            }
            for (std::size_t j = 0; j <= N; ++j) {
                if (used[j]) { u[p[j]] += delta; v[j] -= delta; } else minv[j] -= delta;
            }
            j0 = j1;
        } while (p[j0] != 0);
        do { const auto j1 = way[j0]; p[j0] = p[j1]; j0 = j1; } while (j0);
    }
    VisitMatching out;
    std::vector<bool> firstMatched(m, false), secondMatched(n, false);
    for (std::size_t j = 1; j <= N; ++j) {
        const auto i = p[j] - 1, jj = j - 1;
        if (i < m && jj < n && cost[i * N + jj] <= 1) { out.matches.emplace_back(i, jj); firstMatched[i] = secondMatched[jj] = true; }
    }
    std::sort(out.matches.begin(), out.matches.end());
    for (std::size_t i = 0; i < m; ++i) if (!firstMatched[i]) out.unmatchedFirst.push_back(i);
    for (std::size_t j = 0; j < n; ++j) if (!secondMatched[j]) out.unmatchedSecond.push_back(j);
    return out;
}

}  // namespace dynsyslab::core
