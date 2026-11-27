// Compactly supported radial fields with an optional affine tail (FieldModels.hpp). Ported from
// tools/experiments/miscellaneous_computes/compact_field.py. Weights: Wendland's C^2 function
// (1 - r)^4 (1 + 4 r) of the distance to each of the k nearest centers over the (k+1)-th
// (Wendland 1995, Adv. Comput. Math. 4:389-396), normalized to a partition of unity (Shepard 1968).
#include "FieldModelsInternal.hpp"

#include <algorithm>
#include <cmath>
#include <memory>
#include <set>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

struct Compact {
    std::vector<double> centers, coefficients, linear;  // C x n, C x m, n x m
    std::size_t n{}, m{}, k{};
    double scale{};
    [[nodiscard]] std::vector<double> transformed(const double* x) const {
        std::vector<double> z(x, x + n);
        if (scale > 0) for (auto& v : z) v = std::asinh(v / scale);
        return z;
    }
    // (center index, weight) of the k nearest centers.
    [[nodiscard]] std::vector<std::pair<std::size_t, double>> features(const double* x) const {
        const auto z = transformed(x);
        const auto hits = detail::nearest(centers, n, z.data(), k + 1);
        const auto w = detail::wendlandWeights(hits);
        std::vector<std::pair<std::size_t, double>> out;
        for (std::size_t i = 0; i < w.size(); ++i) out.emplace_back(hits[i].second, w[i]);
        return out;
    }
    [[nodiscard]] std::vector<double> operator()(const double* x) const {
        std::vector<double> y(m, 0);
        for (const auto& [c, w] : features(x)) for (std::size_t j = 0; j < m; ++j) y[j] += w * coefficients[c * m + j];
        for (std::size_t i = 0; i < n; ++i) for (std::size_t j = 0; j < m; ++j) y[j] += x[i] * linear[i * m + j];
        return y;
    }
};

Compact read(const json::Value& d) {
    Compact c;
    c.n = std::size_t(d.at("dimension").number()); c.m = std::size_t(d.at("outputs").number());
    c.k = std::size_t(d.at("neighbors").number()); c.scale = d.at("asinh_scale").number(0);
    c.centers = detail::jsonNumbers(d.at("centers")); c.coefficients = detail::jsonNumbers(d.at("coefficients"));
    c.linear = detail::jsonNumbers(d.at("linear"));
    const auto C = c.n ? c.centers.size() / c.n : 0;
    if (!c.n || !c.m || !c.k || c.k >= C || c.coefficients.size() != C * c.m || c.linear.size() != c.n * c.m || c.scale < 0)
        throw std::invalid_argument("compact_radial: inconsistent definition");
    return c;
}

}  // namespace

namespace detail {
Observable compactRadialField(const json::Value& d) {
    const auto c = std::make_shared<const Compact>(read(d));
    return Observable({"compact_radial", "compact radial field", "", "model", {}, false}, c->n, c->m,
                      [c](const double* x) { const auto y = (*c)(x); return Observable::Values(y.begin(), y.end()); });
}
}  // namespace detail

json::Value fitCompactRadialField(const std::vector<FieldRecord>& records, const CompactRadialOptions& o) {
    if (records.empty() || !o.stride || !o.neighbors || o.asinhScale < 0 || o.speedScale < 0)
        throw std::invalid_argument("compact_radial fit: records, stride, neighbors and scales");
    Compact c;
    c.n = records[0].dimension; c.m = records[0].outputs; c.k = o.neighbors; c.scale = o.asinhScale;
    std::set<std::vector<double>> unique;
    for (const auto& r : records) {
        if (r.dimension != c.n || r.outputs != c.m) throw std::invalid_argument("compact_radial fit: record dimensions disagree");
        for (std::size_t k = 0; k < r.samples(); k += o.stride) unique.emplace(r.state.begin() + std::ptrdiff_t(k * c.n), r.state.begin() + std::ptrdiff_t((k + 1) * c.n));
    }
    std::vector<double> raw;
    for (const auto& row : unique) raw.insert(raw.end(), row.begin(), row.end());
    const auto C = unique.size();
    if (c.k >= C) throw std::invalid_argument("compact_radial fit: fewer centers than neighbours");
    for (std::size_t i = 0; i < C; ++i) { const auto z = c.transformed(raw.data() + i * c.n); c.centers.insert(c.centers.end(), z.begin(), z.end()); }
    const auto tail = o.affineTail ? c.n : 0, F = C + tail, N = F + c.m;
    std::vector<double> gram(N * N, 0);
    for (const auto& r : records) {
        const auto quadrature = o.timeWeighted ? detail::trapezoidWeights(r.time) : std::vector<double>();
        for (std::size_t k = 0; k < r.samples(); ++k) {
            const double* x = r.state.data() + k * c.n;
            const double* y = r.label.data() + k * c.m;
            double weight = o.timeWeighted ? std::sqrt(quadrature[k]) : 1.0;
            if (o.speedScale > 0) {
                double speed = 0; for (std::size_t j = 0; j < c.m; ++j) speed += y[j] * y[j];
                weight /= std::hypot(1.0, std::sqrt(speed) / o.speedScale);
            }
            std::vector<std::pair<std::size_t, double>> row;  // sparse row of [features, tail, targets]
            const auto features = c.features(x);
            for (const auto& [i, w] : features) row.emplace_back(i, w * weight);
            if (tail) {
                std::vector<double> mean(c.n, 0);
                for (const auto& [i, w] : features) for (std::size_t j = 0; j < c.n; ++j) mean[j] += w * raw[i * c.n + j];
                for (std::size_t j = 0; j < c.n; ++j) row.emplace_back(C + j, (x[j] - mean[j]) * weight);
            }
            for (std::size_t j = 0; j < c.m; ++j) row.emplace_back(F + j, y[j] * weight);
            for (const auto& [a, va] : row) for (const auto& [b, vb] : row) gram[a * N + b] += va * vb;
        }
    }
    const auto coef = thresholdedLeastSquares(gram, F, c.m, o.ridge, o.threshold);  // F x m
    c.linear.assign(c.n * c.m, 0);
    if (tail) for (std::size_t i = 0; i < c.n; ++i) for (std::size_t j = 0; j < c.m; ++j) c.linear[i * c.m + j] = coef[(C + i) * c.m + j];
    c.coefficients.assign(C * c.m, 0);
    for (std::size_t a = 0; a < C; ++a)
        for (std::size_t j = 0; j < c.m; ++j) {
            double v = coef[a * c.m + j];
            for (std::size_t i = 0; i < c.n; ++i) v -= raw[a * c.n + i] * c.linear[i * c.m + j];
            c.coefficients[a * c.m + j] = v;
        }
    return json::Value::Object{{"kind", "compact_radial"}, {"dimension", double(c.n)}, {"outputs", double(c.m)},
        {"neighbors", double(c.k)}, {"asinh_scale", c.scale}, {"affine_tail", o.affineTail},
        {"centers", detail::numbersJson(c.centers)}, {"coefficients", detail::numbersJson(c.coefficients)},
        {"linear", detail::numbersJson(c.linear)}};
}

}  // namespace dynsyslab::core
