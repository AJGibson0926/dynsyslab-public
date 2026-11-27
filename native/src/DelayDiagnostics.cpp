// Delay pairs, delay consistency, locked basis extension and the endpoint decoder
// (core/DelayDiagnostics.hpp, which cites the methods).
#include <dynsyslab/core/DelayDiagnostics.hpp>
#include <dynsyslab/numerics/MatrixFunctions.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace dynsyslab::core {

DelayPairs scalarDelayPairs(const std::vector<double>& p, const std::vector<std::size_t>& ends,
                            const std::vector<std::size_t>& lags, double dt, std::size_t h) {
    if (!(dt > 0) || !std::isfinite(dt) || h < 1 || ends.empty() || lags.empty())
        throw std::invalid_argument("delay pairs need ends, lags, a positive step and half-width");
    const auto longest = *std::max_element(lags.begin(), lags.end());
    for (const auto e : ends)
        if (e >= p.size() || e < longest + 2 * h) throw std::invalid_argument("each delay pair must lie inside its record");
    DelayPairs out;
    out.ends = ends.size(); out.lags = lags.size(); out.latency = double(h) * dt;
    out.raw.resize(ends.size() * 2 * lags.size()); out.state.resize(out.raw.size()); out.curvature.resize(ends.size() * lags.size());
    const double delta = double(h) * dt;
    for (std::size_t r = 0; r < ends.size(); ++r)
        for (std::size_t l = 0; l < lags.size(); ++l) {
            const auto latest = ends[r] - lags[l];
            const double older = p[latest - 2 * h], newer = p[latest], middle = p[latest - h];
            const auto at = r * 2 * lags.size() + 2 * l;
            out.raw[at] = older; out.raw[at + 1] = newer;
            out.state[at] = 0.5 * (newer + older); out.state[at + 1] = (newer - older) / (2 * delta);
            out.curvature[r * lags.size() + l] = (newer - 2 * middle + older) / (delta * delta);
        }
    return out;
}

namespace {
// Cubic Hermite value and derivative at tau from (t, y, m) on an increasing clock.
std::pair<double, double> hermite(const std::vector<double>& t, const std::vector<double>& x, std::size_t width, double tau) {
    auto k = std::size_t(std::upper_bound(t.begin(), t.end(), tau) - t.begin());
    k = k == 0 ? 0 : std::min(k - 1, t.size() - 2);
    const double h = t[k + 1] - t[k], s = (tau - t[k]) / h;
    const double y0 = x[k * width], m0 = x[k * width + 1], y1 = x[(k + 1) * width], m1 = x[(k + 1) * width + 1];
    const double s2 = s * s, s3 = s2 * s;
    const double value = (2 * s3 - 3 * s2 + 1) * y0 + (s3 - 2 * s2 + s) * h * m0 + (-2 * s3 + 3 * s2) * y1 + (s3 - s2) * h * m1;
    const double rate = ((6 * s2 - 6 * s) * y0 + (-6 * s2 + 6 * s) * y1) / h + (3 * s2 - 4 * s + 1) * m0 + (3 * s2 - 2 * s) * m1;
    return {value, rate};
}
}  // namespace

DelayConsistency delayConsistency(const std::vector<double>& t, const std::vector<double>& x,
                                  const std::vector<double>& r, const std::vector<double>& lags, std::size_t stride) {
    const auto n = t.size(), width = 2 * lags.size();
    bool valid = n >= 2 && !lags.empty() && lags[0] == 0.0 && stride >= 1 && x.size() == n * width && r.size() == x.size();
    for (std::size_t k = 1; valid && k < n; ++k) valid = t[k] > t[k - 1];
    for (const auto l : lags) valid = valid && l >= 0 && std::isfinite(l);
    for (const auto v : x) valid = valid && std::isfinite(v);
    for (const auto v : r) valid = valid && std::isfinite(v);
    if (!valid) throw std::invalid_argument("delay consistency needs an increasing clock, finite aligned (value, rate) taps and lags from 0");
    const double nan = std::numeric_limits<double>::quiet_NaN();
    DelayConsistency out;
    for (std::size_t j = 1; j < lags.size(); ++j) {
        out.lags.push_back(lags[j]);
        std::vector<std::size_t> ids;
        for (std::size_t k = 0; k < n; ++k) if (t[k] - lags[j] >= t[0]) ids.push_back(k);
        std::vector<std::size_t> kept;
        for (std::size_t i = 0; i < ids.size(); i += stride) kept.push_back(ids[i]);
        if (kept.size() < 2) { out.value.push_back(nan); out.rate.push_back(nan); continue; }
        double e[2] = {0, 0}, d[2] = {0, 0};
        for (const auto k : kept) {
            const auto [value, rate] = hermite(t, x, width, t[k] - lags[j]);
            const double dv = value - x[k * width + 2 * j], dr = rate - x[k * width + 2 * j + 1];
            e[0] += dv * dv; e[1] += dr * dr;
            d[0] += r[k * width + 2 * j] * r[k * width + 2 * j]; d[1] += r[k * width + 2 * j + 1] * r[k * width + 2 * j + 1];
        }
        out.value.push_back(d[0] > 0 ? std::sqrt(e[0] / d[0]) : nan);
        out.rate.push_back(d[1] > 0 ? std::sqrt(e[1] / d[1]) : nan);
    }
    return out;
}

LockedBasis extendLockedBasis(const std::vector<double>& G, std::size_t n, const std::vector<double>& U0,
                              const std::vector<double>& s0, std::size_t rank) {
    const auto p = s0.size();
    if (G.size() != n * n || U0.size() != n * p || p > rank || rank > n || rank == 0)
        throw std::invalid_argument("locked basis needs an n x n Gram matrix, an n x p prefix and p <= rank <= n");
    for (std::size_t a = 0; a < p; ++a) {
        if (!(s0[a] > 0)) throw std::invalid_argument("locked singular values must be positive");
        double residual = 0;
        for (std::size_t i = 0; i < n; ++i) {
            double gu = 0;
            for (std::size_t k = 0; k < n; ++k) gu += G[i * n + k] * U0[k * p + a];
            residual += std::pow(gu - s0[a] * s0[a] * U0[i * p + a], 2);
        }
        for (std::size_t b = 0; b < p; ++b) {
            double dot = 0;
            for (std::size_t i = 0; i < n; ++i) dot += U0[i * p + a] * U0[i * p + b];
            if (std::abs(dot - (a == b ? 1.0 : 0.0)) > 1e-8) throw std::invalid_argument("locked prefix is not orthonormal");
        }
        if (std::sqrt(residual) / (s0[a] * s0[a]) > 1e-6)
            throw std::invalid_argument("frozen eigenpairs do not match this Gram matrix");
    }
    LockedBasis out;
    out.modes.assign(n * rank, 0.0);
    out.singularValues = s0;
    for (std::size_t i = 0; i < n; ++i) for (std::size_t a = 0; a < p; ++a) out.modes[i * rank + a] = U0[i * p + a];
    if (rank > p) {
        std::vector<double> P(n * n, 0.0), PG(n * n, 0.0), PGP(n * n, 0.0);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t k = 0; k < n; ++k) {
                double v = i == k ? 1.0 : 0.0;
                for (std::size_t a = 0; a < p; ++a) v -= U0[i * p + a] * U0[k * p + a];
                P[i * n + k] = v;
            }
        for (std::size_t i = 0; i < n; ++i) for (std::size_t m = 0; m < n; ++m) {
            const double pim = P[i * n + m];
            if (pim != 0.0) for (std::size_t k = 0; k < n; ++k) PG[i * n + k] += pim * G[m * n + k];
        }
        for (std::size_t i = 0; i < n; ++i) for (std::size_t m = 0; m < n; ++m) {
            const double v = PG[i * n + m];
            if (v != 0.0) for (std::size_t k = 0; k < n; ++k) PGP[i * n + k] += v * P[m * n + k];
        }
        for (std::size_t i = 0; i < n; ++i) for (std::size_t k = 0; k < i; ++k)
            PGP[i * n + k] = PGP[k * n + i] = 0.5 * (PGP[i * n + k] + PGP[k * n + i]);
        const auto eig = numerics::symmetricEigensystem(PGP, n);  // ascending
        const double top = std::max(eig.eigenvalues.back(), p ? s0[0] * s0[0] : 0.0);
        for (std::size_t c = 0; c < rank - p; ++c) {
            const auto column = n - 1 - c;
            const double lambda = eig.eigenvalues[column];
            if (lambda <= std::max(1e-14, std::numeric_limits<double>::epsilon() * double(n)) * top)
                throw std::invalid_argument("requested rank exceeds the numerical rank");
            std::size_t pivot = 0;
            for (std::size_t i = 1; i < n; ++i)
                if (std::abs(eig.eigenvectors[i * n + column]) > std::abs(eig.eigenvectors[pivot * n + column])) pivot = i;
            const double sign = eig.eigenvectors[pivot * n + column] < 0 ? -1.0 : 1.0;
            for (std::size_t i = 0; i < n; ++i) out.modes[i * rank + p + c] = sign * eig.eigenvectors[i * n + column];
            out.singularValues.push_back(std::sqrt(lambda));
        }
    }
    const auto& U = out.modes; const auto& s = out.singularValues;
    for (std::size_t a = 0; a < rank; ++a)
        for (std::size_t b = 0; b < rank; ++b) {
            double v = 0;
            for (std::size_t i = 0; i < n; ++i) {
                double gu = 0;
                for (std::size_t k = 0; k < n; ++k) gu += G[i * n + k] * U[k * rank + b];
                v += U[i * rank + a] * gu;
            }
            out.gramError = std::max(out.gramError, std::abs(v / (s[a] * s[b]) - (a == b ? 1.0 : 0.0)));
        }
    return out;
}

std::vector<double> historyEndpointDecoder(const HankelLeftBasis& basis, std::size_t rank, std::size_t channels,
                                           std::size_t channel, double gain) {
    if (!rank || rank > basis.rank || !channels || channel >= channels || basis.rows % channels
        || basis.leftSingularVectors.size() != basis.rows * basis.rank || !(gain > 0))
        throw std::invalid_argument("invalid endpoint decoder rank, channel or gain");
    const auto row = basis.rows - channels + channel;
    std::vector<double> d(rank);
    for (std::size_t i = 0; i < rank; ++i) d[i] = basis.leftSingularVectors[row * basis.rank + i] * basis.singularValues[i] / gain;
    return d;
}

}  // namespace dynsyslab::core
