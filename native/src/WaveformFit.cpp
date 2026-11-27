// Filtered-library waveform fit and the waveform objective (FieldModels.hpp). Ported from
// tools/experiments/miscellaneous_computes/waveform_fit.py and waveform_objective.py.
// The critically damped filter x'' + 2 g x' + g^2 x = g^2 phi is discretized exactly for inputs
// held at interval midpoints (zero-order-hold discretization; Franklin, Powell and Workman 1998,
// "Digital Control of Dynamic Systems", chapter 6), so the regression is linear in the library
// coefficients (output-error identification on filtered regressors; Ljung 1999, section 7.6).
#include "FieldModelsInternal.hpp"

#include <dynsyslab/core/ModelLibraries.hpp>
#include <dynsyslab/core/SindyDynamics.hpp>
#include <dynsyslab/core/TimeSeries.hpp>
#include <dynsyslab/numerics/Svd.hpp>

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

// Position and velocity responses (samples x columns each) of every column of `features`.
std::pair<std::vector<double>, std::vector<double>> filteredLibrary(const std::vector<double>& t, const std::vector<double>& features,
                                                                    std::size_t columns, double gamma) {
    const auto S = t.size();
    if (S < 2 || features.size() != S * columns || !std::isfinite(gamma) || gamma <= 0)
        throw std::invalid_argument("damped drive: a uniform clock, aligned features and a positive damping rate");
    const double h = (t.back() - t.front()) / double(S - 1);
    for (std::size_t k = 1; k < S; ++k)
        if (!(h > 0) || std::abs(t[k] - t[k - 1] - h) > 1e-6 * h) throw std::invalid_argument("damped drive: the clock must be uniform");
    const double x = gamma * h, q = std::exp(-x);
    const double h0 = x > 1e-4 ? (-std::expm1(-x) - x * q) : x * x * (0.5 - x / 3 + x * x / 8), h1 = gamma * x * q;
    const double A00 = q * (1 + x), A01 = h * q, A10 = -gamma * x * q, A11 = (1 - x) * q;
    const double num[2][2] = {{h0, A01 * h1 - A11 * h0}, {h1, A10 * h0 - A00 * h1}};
    std::vector<double> out[2] = {std::vector<double>(S * columns, 0), std::vector<double>(S * columns, 0)};
    for (int which = 0; which < 2; ++which)
        for (std::size_t c = 0; c < columns; ++c) {
            double y1 = 0, y2 = 0, m1 = 0;  // y[k-1], y[k-2], mid[k-1]
            for (std::size_t k = 0; k + 1 < S; ++k) {
                const double mid = 0.5 * (features[k * columns + c] + features[(k + 1) * columns + c]);
                const double y = num[which][0] * mid + (k ? num[which][1] * m1 : 0) + 2 * q * y1 - q * q * y2;
                out[which][(k + 1) * columns + c] = y;
                y2 = y1; y1 = y; m1 = mid;
            }
        }
    return {out[0], out[1]};
}

}  // namespace

json::Value fitDampedDrive(const std::vector<DriveRecord>& records, double gamma, const MapFitOptions& o) {
    if (records.empty()) throw std::invalid_argument("damped drive: no records");
    const auto m = records[0].inputs, p = records[0].outputs;
    std::vector<double> all, index;
    for (const auto& r : records) {
        const auto S = r.time.size();
        if (r.inputs != m || r.outputs != p || r.input.size() != S * m || r.position.size() != S * p || r.velocity.size() != S * p)
            throw std::invalid_argument("damped drive: aligned inputs and position/velocity labels");
        all.insert(all.end(), r.input.begin(), r.input.end());
        for (std::size_t k = 0; k < S; ++k) index.push_back(double(index.size()));
    }
    ModelLibraryOptions lo;
    lo.degree = o.degree; lo.centers = o.centers; lo.harmonics = o.harmonics; lo.rationalWidths = o.rationalWidths;
    const auto library = modelLibrary(o.family, TimeSeries(index, all, m), lo);
    const auto F = library.outputDimension(), N = F + p;
    std::vector<double> gram(N * N, 0);
    for (const auto& r : records) {
        const auto S = r.time.size();
        std::vector<double> theta;
        theta.reserve(S * F);
        for (std::size_t k = 0; k < S; ++k) { const auto row = library.evaluateReal(r.input.data() + k * m); theta.insert(theta.end(), row.begin(), row.end()); }
        const auto [pos, vel] = filteredLibrary(r.time, theta, F, gamma);
        for (int which = 0; which < 2; ++which) {
            const auto& design = which ? vel : pos;
            const auto& observed = which ? r.velocity : r.position;
            double energy = 0;
            for (const auto v : observed) energy += v * v;
            const double scale = energy > 0 ? std::sqrt(energy / double(observed.size())) : 1.0;
            std::vector<double> row(N);
            for (std::size_t k = 0; k < S; ++k) {
                const double dt = r.time[k] - r.time[0], e = std::exp(-gamma * dt);
                for (std::size_t s = 0; s < F; ++s) row[s] = design[k * F + s] / scale;
                for (std::size_t j = 0; j < p; ++j) {
                    const double p0 = r.position[j], v0 = r.velocity[j];
                    const double transient = which ? e * ((1 - gamma * dt) * v0 - gamma * gamma * dt * p0) : e * ((1 + gamma * dt) * p0 + dt * v0);
                    row[F + j] = (observed[k * p + j] - transient) / scale;
                }
                for (std::size_t a = 0; a < N; ++a) for (std::size_t b = 0; b < N; ++b) gram[a * N + b] += row[a] * row[b];
            }
        }
    }
    auto coef = thresholdedLeastSquares(gram, F, p, o.ridge, o.threshold);
    for (auto& c : coef) c *= gamma * gamma;
    const SindyDynamicsModel model(library, coef, p);
    return json::Value::Object{{"kind", "polynomial_model"}, {"model", model.definition()}, {"family", o.family}};
}

std::vector<double> waveformWeights(const std::vector<double>& t, const std::vector<double>& y, std::size_t channels,
                                    const std::vector<double>& globalWeights, const std::vector<double>& localWeights, double duration,
                                    const std::vector<double>& floors, std::vector<double> globalFloors) {
    const auto S = t.size();
    if (!channels || y.size() != S * channels || S < 2 || !std::isfinite(duration) || duration <= 0)
        throw std::invalid_argument("waveform weights: aligned finite channels and a positive interval");
    const auto broadcast = [channels](const std::vector<double>& v) {
        if (v.size() == 1) return std::vector<double>(channels, v[0]);
        if (v.size() != channels) throw std::invalid_argument("waveform weights: one value or one per channel");
        return v;
    };
    const auto gw = broadcast(globalWeights), lw = broadcast(localWeights), floor = broadcast(floors);
    const auto gf = globalFloors.empty() ? floor : broadcast(globalFloors);
    for (std::size_t c = 0; c < channels; ++c)
        if (!(gw[c] >= 0) || !(lw[c] >= 0) || !(floor[c] > 0) || !(gf[c] > 0) || !std::isfinite(gw[c] + lw[c] + floor[c] + gf[c]))
            throw std::invalid_argument("waveform weights: nonnegative weights and strictly positive RMS floors");
    for (const auto v : y) if (!std::isfinite(v)) throw std::invalid_argument("waveform weights: finite observations");
    auto mass = detail::trapezoidWeights(t);
    double total = 0;
    for (const auto m : mass) total += m;
    for (auto& m : mass) m /= total;
    std::vector<double> global(channels, 0);
    for (std::size_t k = 0; k < S; ++k) for (std::size_t c = 0; c < channels; ++c) global[c] += mass[k] * y[k * channels + c] * y[k * channels + c];
    for (std::size_t c = 0; c < channels; ++c) global[c] = std::max(global[c], gf[c] * gf[c]);
    std::map<long long, std::vector<std::size_t>> blocks;
    for (std::size_t k = 0; k < S; ++k) blocks[(long long)std::floor((t[k] - t[0]) / duration)].push_back(k);
    std::vector<double> local(S * channels);
    for (const auto& [block, rows] : blocks) {
        double w = 0;
        std::vector<double> sum(channels, 0);
        for (const auto k : rows) { w += mass[k]; for (std::size_t c = 0; c < channels; ++c) sum[c] += mass[k] * y[k * channels + c] * y[k * channels + c]; }
        for (std::size_t c = 0; c < channels; ++c) {
            const double value = std::max(w > 0 ? sum[c] / w : 0.0, floor[c] * floor[c]);
            for (const auto k : rows) local[k * channels + c] = value;
        }
    }
    std::vector<double> out(S * channels);
    for (std::size_t k = 0; k < S; ++k)
        for (std::size_t c = 0; c < channels; ++c)
            out[k * channels + c] = std::sqrt(mass[k] * (gw[c] / global[c] + lw[c] / local[k * channels + c]));
    return out;
}

std::vector<double> resolvedDirections(const std::vector<double>& x, std::size_t rows, const std::vector<double>& v, std::size_t coefficients,
                                       std::size_t p, std::vector<double>& responses, double relativeCutoff, double singularFloor, double scale) {
    if (!rows || !p || !coefficients || x.size() != rows * p || v.size() != coefficients * p || !(relativeCutoff >= 0 && relativeCutoff < 1)
        || !(singularFloor > 0) || !(scale > 0) || !std::isfinite(singularFloor + scale))
        throw std::invalid_argument("resolved directions: invalid training response or coefficient subspace");
    for (const auto e : x) if (!std::isfinite(e)) throw std::invalid_argument("resolved directions: finite response");
    for (const auto e : v) if (!std::isfinite(e)) throw std::invalid_argument("resolved directions: finite directions");
    std::vector<double> norms(p, 0);
    for (std::size_t r = 0; r < rows; ++r) for (std::size_t j = 0; j < p; ++j) norms[j] += x[r * p + j] * x[r * p + j];
    double largest = 0;
    for (auto& n : norms) { n = std::sqrt(n); largest = std::max(largest, n); }
    std::vector<std::size_t> active;
    for (std::size_t j = 0; j < p; ++j) if (norms[j] > largest * 1e-10) active.push_back(j);
    if (active.empty()) throw std::invalid_argument("resolved directions: training data resolve no supplied direction");
    const auto a = active.size();
    std::vector<double> scaled(rows * a);
    for (std::size_t r = 0; r < rows; ++r) for (std::size_t i = 0; i < a; ++i) scaled[r * a + i] = x[r * p + active[i]] / norms[active[i]];
    const auto svd = numerics::singularValueDecomposition(scaled, rows, a);
    std::vector<std::size_t> keep;
    for (std::size_t k = 0; k < svd.rank; ++k) if (svd.singularValues[k] > svd.singularValues[0] * relativeCutoff) keep.push_back(k);
    const auto K = keep.size();
    std::vector<double> w(a * K);  // active x kept
    for (std::size_t i = 0; i < a; ++i)
        for (std::size_t q = 0; q < K; ++q)
            w[i * K + q] = svd.rightTranspose[keep[q] * a + i] / std::max(svd.singularValues[keep[q]], singularFloor) * scale / norms[active[i]];
    std::vector<double> directions(coefficients * K, 0);
    responses.assign(rows * K, 0);
    for (std::size_t c = 0; c < coefficients; ++c)
        for (std::size_t i = 0; i < a; ++i) for (std::size_t q = 0; q < K; ++q) directions[c * K + q] += v[c * p + active[i]] * w[i * K + q];
    for (std::size_t r = 0; r < rows; ++r)
        for (std::size_t i = 0; i < a; ++i) for (std::size_t q = 0; q < K; ++q) responses[r * K + q] += x[r * p + active[i]] * w[i * K + q];
    return directions;
}

}  // namespace dynsyslab::core
