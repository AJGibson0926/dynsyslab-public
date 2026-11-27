// Neighbour readouts and oscillation-interval errors (NeighborReadout.hpp). Ported from the
// receiver-reconstruction side-car (lookup.cpp, rings.cpp, signals.cpp maxima).
#include <dynsyslab/core/NeighborReadout.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace dynsyslab::core {
namespace {
bool closer(double da, std::size_t ra, double db, std::size_t rb) { return da == db ? ra < rb : da < db; }
}  // namespace

NeighborReadout::NeighborReadout(std::vector<double> features, std::vector<double> states, std::vector<double> scales,
                                 std::size_t stateDimension)
    : points_(std::move(features)), states_(std::move(states)), scales_(std::move(scales)), stateDimension_(stateDimension) {
    const auto d = scales_.size();
    if (!d || !stateDimension_ || points_.empty() || points_.size() % d || states_.size() % stateDimension_
        || states_.size() / stateDimension_ != points_.size() / d)
        throw std::invalid_argument("neighbour readout: training feature and state shapes disagree");
    for (const double s : scales_) if (!std::isfinite(s) || s <= 0) throw std::invalid_argument("neighbour readout: scales must be positive");
    for (const double s : states_) if (!std::isfinite(s)) throw std::invalid_argument("neighbour readout: nonfinite training state");
    for (std::size_t i = 0; i < points_.size(); ++i) {
        points_[i] = std::asinh(points_[i] / scales_[i % d]);
        if (!std::isfinite(points_[i])) throw std::invalid_argument("neighbour readout: nonfinite training feature");
    }
    order_.resize(points_.size() / d);
    std::iota(order_.begin(), order_.end(), std::size_t(0));
    build(0, order_.size(), 0);
}

void NeighborReadout::build(std::size_t first, std::size_t last, std::size_t axis) {
    if (last - first < 2) return;
    const auto mid = first + (last - first) / 2, d = scales_.size();
    std::nth_element(order_.begin() + std::ptrdiff_t(first), order_.begin() + std::ptrdiff_t(mid), order_.begin() + std::ptrdiff_t(last),
                     [&](auto a, auto b) { const double x = points_[a * d + axis], y = points_[b * d + axis]; return x == y ? a < b : x < y; });
    build(first, mid, (axis + 1) % d);
    build(mid + 1, last, (axis + 1) % d);
}

void NeighborReadout::search(std::size_t first, std::size_t last, std::size_t axis, const std::vector<double>& query, std::size_t k,
                             std::vector<Hit>& hits) const {
    if (first == last) return;
    const auto mid = first + (last - first) / 2, row = order_[mid], d = scales_.size();
    const double* point = points_.data() + row * d;
    double square = 0;
    for (std::size_t j = 0; j < d; ++j) square += (query[j] - point[j]) * (query[j] - point[j]);
    if (hits.size() < k || closer(square, row, hits.back().squared, hits.back().row)) {
        if (hits.size() == k) hits.pop_back();
        const Hit hit{square, row};
        hits.insert(std::lower_bound(hits.begin(), hits.end(), hit,
                                     [](const Hit& a, const Hit& b) { return closer(a.squared, a.row, b.squared, b.row); }), hit);
    }
    const double delta = query[axis] - point[axis];
    const auto next = (axis + 1) % d;
    const auto nearSide = [&](bool low) { if (low) search(first, mid, next, query, k, hits); else search(mid + 1, last, next, query, k, hits); };
    nearSide(delta <= 0);
    if (hits.size() < k || delta * delta <= hits.back().squared) nearSide(!(delta <= 0));
}

NeighborReadout::Estimate NeighborReadout::predict(const double* features, std::size_t k) const {
    if (!features || !k || k > order_.size()) throw std::invalid_argument("neighbour readout: invalid query or neighbour count");
    std::vector<double> query(scales_.size());
    for (std::size_t j = 0; j < query.size(); ++j) {
        query[j] = std::asinh(features[j] / scales_[j]);
        if (!std::isfinite(query[j])) throw std::invalid_argument("neighbour readout: nonfinite query feature");
    }
    std::vector<Hit> hits;
    hits.reserve(k);
    search(0, order_.size(), 0, query, k, hits);
    Estimate out;
    const double nearest = std::max(1e-24, hits.front().squared);
    out.nearestDistance = std::sqrt(hits.front().squared);
    double sum = 0;
    for (const auto& h : hits) {
        out.rows.push_back(h.row);
        out.weights.push_back(nearest / std::max(1e-24, h.squared));
        sum += out.weights.back();
    }
    out.values.assign(stateDimension_, 0);
    for (std::size_t i = 0; i < hits.size(); ++i) {
        out.weights[i] /= sum;
        for (std::size_t j = 0; j < stateDimension_; ++j) out.values[j] += out.weights[i] * states_[hits[i].row * stateDimension_ + j];
    }
    return out;
}

std::vector<double> oscillationMaxima(const std::vector<double>& t, const std::vector<double>& x, std::size_t width, std::size_t velocity) {
    if (!width || velocity >= width || x.size() != t.size() * width || t.size() < 2)
        throw std::invalid_argument("oscillation maxima: shapes");
    std::vector<double> out;
    if (x[velocity] == 0 && x[width + velocity] < 0) out.push_back(t.front());  // released from rest
    for (std::size_t i = 1; i < t.size(); ++i) {
        const double a = x[(i - 1) * width + velocity], b = x[i * width + velocity];
        if (a > 0 && b <= 0) out.push_back(t[i - 1] + a / (a - b) * (t[i] - t[i - 1]));
    }
    return out;
}

IntervalErrors intervalErrors(const std::vector<double>& t, const std::vector<double>& x, const std::vector<double>& y, std::size_t d,
                              const std::vector<double>& bounds, std::vector<double> center) {
    if (!d || t.empty() || x.size() != t.size() * d || y.size() != x.size() || bounds.size() < 2)
        throw std::invalid_argument("interval errors need aligned reference and prediction and at least two bounds");
    const auto increasing = [](const std::vector<double>& v) {
        return std::all_of(v.begin(), v.end(), [](double z) { return std::isfinite(z); })
            && std::adjacent_find(v.begin(), v.end(), std::greater_equal<double>()) == v.end();
    };
    if (!increasing(t) || !increasing(bounds)) throw std::invalid_argument("interval errors: clocks and bounds must increase");
    if (center.empty()) center.assign(d, 0);
    if (center.size() != d) throw std::invalid_argument("interval errors: center size");
    const auto n = bounds.size() - 1;
    IntervalErrors out;
    out.width = d;
    out.count.assign(n, 0); out.partial.assign(n, 0);
    out.errorSquares.assign(n * d, 0); out.referenceSquares.assign(n * d, 0); out.relativeRmsPercent.assign(n * d, 0);
    for (std::size_t i = 0; i < n; ++i) {
        const auto first = std::size_t(std::lower_bound(t.begin(), t.end(), bounds[i]) - t.begin());
        const auto end = std::size_t(std::lower_bound(t.begin(), t.end(), bounds[i + 1]) - t.begin());
        out.count[i] = double(end - first);
        out.partial[i] = bounds[i] < t.front() || bounds[i + 1] > t.back();
        for (auto k = first; k < end; ++k)
            for (std::size_t j = 0; j < d; ++j) {
                const double e = y[k * d + j] - x[k * d + j], r = x[k * d + j] - center[j];
                if (!std::isfinite(e) || !std::isfinite(r)) throw std::invalid_argument("interval errors: nonfinite values");
                out.errorSquares[i * d + j] += e * e;
                out.referenceSquares[i * d + j] += r * r;
            }
        for (std::size_t j = 0; j < d; ++j) {
            const auto at = i * d + j;
            out.relativeRmsPercent[at] = out.referenceSquares[at] > 0 ? 100 * std::sqrt(out.errorSquares[at] / out.referenceSquares[at])
                                                                     : std::numeric_limits<double>::quiet_NaN();
        }
    }
    return out;
}

}  // namespace dynsyslab::core
