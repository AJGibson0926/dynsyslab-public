// Local tangent-chart quadratic fields (FieldModels.hpp). Ported from
// tools/experiments/miscellaneous_computes/local_field.py and its native evaluator local_rhs.cpp.
// Around each training center a Gaussian-weighted local PCA gives a two-dimensional tangent chart
// (local tangent space alignment; Zhang and Zha 2004, SIAM J. Sci. Comput. 26:313-338); a sparse
// quadratic field in chart coordinates is fitted there. Queries blend the nearest centers' fields
// with Wendland weights (Wendland 1995) and damp the component transverse to each chart.
#include "FieldModelsInternal.hpp"

#include <dynsyslab/numerics/MatrixFunctions.hpp>

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

void quadratic(double a, double b, double* f) { f[0] = 1; f[1] = a; f[2] = b; f[3] = a * a; f[4] = a * b; f[5] = b * b; }

struct Local {
    std::size_t n{}, count{}, neighbors{};
    double decay{};
    std::vector<double> center, mean, P, Q, coef, scale;  // count x n, count x n, count x n x 2 (P, Q), count x 6 x n, n
    [[nodiscard]] std::vector<double> operator()(const double* b) const {
        std::vector<double> x(n);
        for (std::size_t i = 0; i < n; ++i) x[i] = b[i] / scale[i];
        const auto hits = detail::nearest(center, n, x.data(), neighbors + 1);
        const auto w = detail::wendlandWeights(hits);
        std::vector<double> out(n, 0), delta(n);
        for (std::size_t h = 0; h < w.size(); ++h) {
            const auto c = hits[h].second;
            for (std::size_t i = 0; i < n; ++i) delta[i] = x[i] - mean[c * n + i];
            double a = 0, bb = 0, ta = 0, tb = 0;
            for (std::size_t i = 0; i < n; ++i) {
                a += delta[i] * P[(c * n + i) * 2]; bb += delta[i] * P[(c * n + i) * 2 + 1];
                ta += delta[i] * Q[(c * n + i) * 2]; tb += delta[i] * Q[(c * n + i) * 2 + 1];
            }
            double f[6];
            quadratic(a, bb, f);
            for (std::size_t i = 0; i < n; ++i) {
                double v = 0;
                for (std::size_t t = 0; t < 6; ++t) v += f[t] * coef[(c * 6 + t) * n + i];
                const double tangent = ta * Q[(c * n + i) * 2] + tb * Q[(c * n + i) * 2 + 1];
                v -= decay * (delta[i] - tangent);
                out[i] += w[h] * v;
            }
        }
        for (std::size_t i = 0; i < n; ++i) out[i] *= scale[i];
        return out;
    }
};

}  // namespace

namespace detail {
Observable localTangentField(const json::Value& d) {
    auto l = std::make_shared<Local>();
    l->n = std::size_t(d.at("dimension").number()); l->neighbors = std::size_t(d.at("neighbors").number()); l->decay = d.at("decay").number();
    l->center = jsonNumbers(d.at("center")); l->mean = jsonNumbers(d.at("mean")); l->P = jsonNumbers(d.at("P"));
    l->Q = jsonNumbers(d.at("Q")); l->coef = jsonNumbers(d.at("coef")); l->scale = jsonNumbers(d.at("scale"));
    l->count = l->n ? l->center.size() / l->n : 0;
    const auto c = l->count, n = l->n;
    if (n < 2 || !l->neighbors || l->neighbors >= c || l->mean.size() != c * n || l->P.size() != c * n * 2 || l->Q.size() != c * n * 2
        || l->coef.size() != c * 6 * n || l->scale.size() != n || l->decay < 0)
        throw std::invalid_argument("local_tangent: inconsistent definition");
    const std::shared_ptr<const Local> field = l;
    return Observable({"local_tangent", "local tangent-chart field", "", "model", {}, false}, n, n,
                      [field](const double* x) { const auto y = (*field)(x); return Observable::Values(y.begin(), y.end()); });
}
}  // namespace detail

json::Value fitLocalTangentField(const std::vector<FieldRecord>& records, const LocalTangentOptions& o) {
    if (records.empty()) throw std::invalid_argument("local_tangent fit: no records");
    const auto n = records[0].dimension;
    std::vector<double> b, v;
    for (const auto& r : records) {
        if (r.dimension != n || r.outputs != n) throw std::invalid_argument("local_tangent fit: records need rates of their own states");
        b.insert(b.end(), r.state.begin(), r.state.end()); v.insert(v.end(), r.label.begin(), r.label.end());
    }
    const auto N = b.size() / n;
    if (n < 2 || !o.stride || o.samples < 6 || o.samples > N || !(o.ridge > 0) || o.threshold < 0)
        throw std::invalid_argument("local_tangent fit: options or too few samples");
    Local l;
    l.n = n; l.neighbors = o.neighbors; l.decay = o.decay; l.scale.assign(n, 0);
    for (std::size_t k = 0; k < N; ++k) for (std::size_t i = 0; i < n; ++i) l.scale[i] += b[k * n + i] * b[k * n + i] / double(N);
    for (auto& s : l.scale) { s = std::sqrt(s); if (s < 1e-14) throw std::invalid_argument("local_tangent fit: a zero-energy coordinate"); }
    std::vector<double> x(N * n), y(N * n);
    for (std::size_t k = 0; k < N; ++k) for (std::size_t i = 0; i < n; ++i) { x[k * n + i] = b[k * n + i] / l.scale[i]; y[k * n + i] = v[k * n + i] / l.scale[i]; }
    for (const auto& r : records)
        for (std::size_t k = 0; k < r.samples(); k += o.stride)
            for (std::size_t i = 0; i < n; ++i) l.center.push_back(r.state[k * n + i] / l.scale[i]);
    l.count = l.center.size() / n;
    for (std::size_t c = 0; c < l.count; ++c) {
        const auto hits = detail::nearest(x, n, l.center.data() + c * n, o.samples);
        const double far = std::max(hits.back().first, 1e-12);
        std::vector<double> w(hits.size()), mean(n, 0);
        double total = 0;
        for (std::size_t h = 0; h < hits.size(); ++h) { w[h] = std::exp(-2.5 * std::pow(hits[h].first / far, 2)); total += w[h]; }
        for (auto& wh : w) wh /= total;
        for (std::size_t h = 0; h < hits.size(); ++h) for (std::size_t i = 0; i < n; ++i) mean[i] += w[h] * x[hits[h].second * n + i];
        std::vector<double> cov(n * n, 0);
        for (std::size_t h = 0; h < hits.size(); ++h)
            for (std::size_t i = 0; i < n; ++i)
                for (std::size_t j = 0; j < n; ++j)
                    cov[i * n + j] += w[h] * (x[hits[h].second * n + i] - mean[i]) * (x[hits[h].second * n + j] - mean[j]);
        const auto eig = numerics::symmetricEigensystem(cov, n);  // ascending
        std::vector<double> P(n * 2), U(n * 2);
        for (std::size_t s = 0; s < 2; ++s) {
            const auto col = n - 2 + s;
            const double norm = std::sqrt(std::max(eig.eigenvalues[col], 1e-7));
            for (std::size_t i = 0; i < n; ++i) { U[i * 2 + s] = eig.eigenvectors[i * n + col]; P[i * 2 + s] = U[i * 2 + s] / norm; }
        }
        // Weighted quadratic regression in chart coordinates, per output with its own support.
        std::vector<double> G(36, 0), rhs(6 * n, 0), energy(n, 0);
        for (std::size_t h = 0; h < hits.size(); ++h) {
            const auto row = hits[h].second;
            double a = 0, bb = 0;
            for (std::size_t i = 0; i < n; ++i) { const double dlt = x[row * n + i] - mean[i]; a += dlt * P[i * 2]; bb += dlt * P[i * 2 + 1]; }
            double f[6];
            quadratic(a, bb, f);
            for (std::size_t s = 0; s < 6; ++s) {
                for (std::size_t t = 0; t < 6; ++t) G[s * 6 + t] += w[h] * f[s] * f[t];
                for (std::size_t i = 0; i < n; ++i) rhs[s * n + i] += f[s] * w[h] * y[row * n + i];
            }
            for (std::size_t i = 0; i < n; ++i) energy[i] += w[h] * y[row * n + i] * y[row * n + i];
        }
        for (std::size_t s = 0; s < 6; ++s) G[s * 6 + s] += o.ridge * std::max(G[s * 6 + s], 1e-10);
        auto coef = detail::choleskySolve(G, 6, rhs, n);
        for (std::size_t i = 0; i < n; ++i) {
            std::vector<std::size_t> keep{0};
            for (std::size_t s = 1; s < 6; ++s)
                if (std::abs(coef[s * n + i]) * std::sqrt(G[s * 6 + s]) > o.threshold * std::sqrt(energy[i])) keep.push_back(s);
            const auto a = keep.size();
            std::vector<double> A(a * a), r(a);
            for (std::size_t p = 0; p < a; ++p) { r[p] = rhs[keep[p] * n + i]; for (std::size_t q = 0; q < a; ++q) A[p * a + q] = G[keep[p] * 6 + keep[q]]; }
            const auto sol = detail::choleskySolve(A, a, r, 1);
            for (std::size_t s = 0; s < 6; ++s) coef[s * n + i] = 0;
            for (std::size_t p = 0; p < a; ++p) coef[keep[p] * n + i] = sol[p];
        }
        l.mean.insert(l.mean.end(), mean.begin(), mean.end());
        l.P.insert(l.P.end(), P.begin(), P.end());
        l.Q.insert(l.Q.end(), U.begin(), U.end());
        l.coef.insert(l.coef.end(), coef.begin(), coef.end());
    }
    if (l.neighbors >= l.count) throw std::invalid_argument("local_tangent fit: fewer centers than neighbours");
    return json::Value::Object{{"kind", "local_tangent"}, {"dimension", double(n)}, {"neighbors", double(l.neighbors)}, {"decay", l.decay},
        {"center", detail::numbersJson(l.center)}, {"mean", detail::numbersJson(l.mean)}, {"P", detail::numbersJson(l.P)},
        {"Q", detail::numbersJson(l.Q)}, {"coef", detail::numbersJson(l.coef)}, {"scale", detail::numbersJson(l.scale)}};
}

}  // namespace dynsyslab::core
