// mode_selection, part 2: score one subset of units.
//
// The subset is realized in scaled real coordinates z = (Re psi, Im psi) / scale per pair
// (psi / scale for a real mode), so A = blockdiag([[s, -w], [w, s]]), B = (Re B, Im B) / scale.
// The readouts C_e and C_x are least-squares fits (with intercepts) on the readout records:
// unforced by default; readout_records selects forced, all, or listed sources.
#include "ModeSelectionInternal.hpp"

#include <dynsyslab/core/BalancedTruncation.hpp>
#include <dynsyslab/core/ControlDesign.hpp>
#include <dynsyslab/core/StateSpace.hpp>
#include <dynsyslab/core/SystemStructure.hpp>
#include <dynsyslab/numerics/MatrixFunctions.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail::mode_selection {
namespace c = dynsyslab::core;
namespace {

std::size_t width(const std::vector<Unit>& units, const std::vector<std::size_t>& subset) {
    std::size_t k = 0; for (auto i : subset) k += units[i].pair ? 2 : 1; return k;
}

void features(const std::vector<Unit>& units, const std::vector<std::size_t>& subset,
              const Record& r, std::size_t i, double* out) {
    std::size_t c0 = 0;
    for (auto u : subset) {
        const auto& unit = units[u];
        const Complex psi = r.psi[i * r.modes + unit.modes[0]] / unit.scale;
        out[c0++] = psi.real();
        if (unit.pair) out[c0++] = psi.imag();
    }
}

// Solve G x = b in place (small, symmetric positive definite after ridge).
std::vector<double> solve(std::vector<double> G, std::vector<double> b, std::size_t n) {
    for (std::size_t k = 0; k < n; ++k) {
        std::size_t p = k;
        for (std::size_t i = k + 1; i < n; ++i) if (std::abs(G[i * n + k]) > std::abs(G[p * n + k])) p = i;
        if (!(std::abs(G[p * n + k]) > 1e-300)) throw std::invalid_argument("readout fit is singular");
        if (p != k) { for (std::size_t j = 0; j < n; ++j) std::swap(G[k * n + j], G[p * n + j]); std::swap(b[k], b[p]); }
        for (std::size_t i = k + 1; i < n; ++i) {
            const double f = G[i * n + k] / G[k * n + k];
            for (std::size_t j = k; j < n; ++j) G[i * n + j] -= f * G[k * n + j];
            b[i] -= f * b[k];
        }
    }
    std::vector<double> x(n);
    for (std::size_t k = n; k-- > 0;) {
        double s = b[k]; for (std::size_t j = k + 1; j < n; ++j) s -= G[k * n + j] * x[j];
        x[k] = s / G[k * n + k];
    }
    return x;
}

double minEigenvalueRatio(const std::vector<double>& M, std::size_t n, double* minimum = nullptr) {
    const auto e = numerics::symmetricEigensystem(M, n).eigenvalues;  // ascending
    if (minimum) *minimum = e.front();
    return e.back() > 0 ? std::max(e.front(), 0.0) / e.back() : 0.0;
}
}  // namespace

Readouts fitReadouts(const std::vector<Unit>& units, const std::vector<std::size_t>& subset,
                     const std::vector<Record>& records, double ridge) {
    const std::size_t k = width(units, subset);
    // Records flagged by readout_records; records built without flags fall back to unforced ones.
    const bool flagged = std::any_of(records.begin(), records.end(), [](const Record& r) { return r.readout; });
    const auto use = [&](const Record& r) { return flagged ? r.readout : !r.forced; };
    std::size_t d = 0, count = 0;
    for (const auto& r : records) if (use(r)) { d = r.stateDimension; count += r.samples; }
    const std::size_t targets = 1 + d;  // echo, then each state component
    std::vector<double> mean(k, 0.0), ymean(targets, 0.0), f(k);
    auto target = [&](const Record& r, std::size_t i, std::size_t t) { return t == 0 ? r.echo[i] : r.states[i * d + t - 1]; };
    for (const auto& r : records) if (use(r))
        for (std::size_t i = 0; i < r.samples; ++i) {
            features(units, subset, r, i, f.data());
            for (std::size_t a = 0; a < k; ++a) mean[a] += f[a];
            for (std::size_t t = 0; t < targets; ++t) ymean[t] += target(r, i, t);
        }
    for (auto& v : mean) v /= double(count);
    for (auto& v : ymean) v /= double(count);
    std::vector<double> G(k * k, 0.0), rhs(k * targets, 0.0), yy(targets, 0.0);
    for (const auto& r : records) if (use(r))
        for (std::size_t i = 0; i < r.samples; ++i) {
            features(units, subset, r, i, f.data());
            for (std::size_t a = 0; a < k; ++a) f[a] -= mean[a];
            for (std::size_t a = 0; a < k; ++a) for (std::size_t b = 0; b < k; ++b) G[a * k + b] += f[a] * f[b];
            for (std::size_t t = 0; t < targets; ++t) {
                const double y = target(r, i, t) - ymean[t];
                yy[t] += y * y;
                for (std::size_t a = 0; a < k; ++a) rhs[t * k + a] += f[a] * y;
            }
        }
    double diag = 0; for (std::size_t a = 0; a < k; ++a) diag = std::max(diag, G[a * k + a]);
    for (std::size_t a = 0; a < k; ++a) G[a * k + a] += ridge * (diag > 0 ? diag : 1.0);
    Readouts out; out.state.assign(d * k, 0.0); out.stateOffset.assign(d, 0.0);
    for (std::size_t t = 0; t < targets; ++t) {
        const auto x = solve(G, {rhs.begin() + std::ptrdiff_t(t * k), rhs.begin() + std::ptrdiff_t((t + 1) * k)}, k);
        double explained = 0, offset = ymean[t];
        for (std::size_t a = 0; a < k; ++a) { explained += x[a] * rhs[t * k + a]; offset -= x[a] * mean[a]; }
        const double residual = std::max(yy[t] - explained, 0.0), sd = std::sqrt(yy[t] / double(count));
        const double error = yy[t] > 0 ? std::sqrt(residual / yy[t]) : 0.0;
        if (t == 0) { out.echo = x; out.echoOffset = offset; out.echoError = error; out.echoStd = sd; }
        else {
            std::copy(x.begin(), x.end(), out.state.begin() + std::ptrdiff_t((t - 1) * k));
            out.stateOffset[t - 1] = offset; out.stateError.push_back(error); out.stateStd.push_back(sd);
        }
    }
    return out;
}

Score scoreSubset(const std::vector<Unit>& units, const std::vector<std::size_t>& subset,
                  const std::vector<Record>& records, const Options& options) {
    Score s;
    try {
        s.readouts = fitReadouts(units, subset, records, options.ridge);
        const std::size_t k = width(units, subset), d = s.readouts.stateStd.size();
        std::vector<double> A(k * k, 0.0), B(k, 0.0), W(k * k, 0.0);
        double dt = 0; for (const auto& r : records) if (!r.forced) { dt = r.dt; break; }
        std::size_t c0 = 0;
        for (auto u : subset) {
            const auto& unit = units[u];
            s.inputMapResidual = std::max(s.inputMapResidual, unit.inputMapResidual);
            const Complex b = unit.inputMap / unit.scale;
            const double q = std::max(unit.processIntensity, 1e-12);
            if (unit.pair) {
                A[c0 * k + c0] = unit.lambda.real(); A[c0 * k + c0 + 1] = -unit.lambda.imag();
                A[(c0 + 1) * k + c0] = unit.lambda.imag(); A[(c0 + 1) * k + c0 + 1] = unit.lambda.real();
                B[c0] = b.real(); B[c0 + 1] = b.imag();
                W[c0 * k + c0] = W[(c0 + 1) * k + c0 + 1] = 0.5 * q;
                c0 += 2;
            } else { A[c0 * k + c0] = unit.lambda.real(); B[c0] = b.real(); W[c0 * k + c0] = q; ++c0; }
        }
        const auto echo = c::makeStateSpace(A, B, s.readouts.echo, {}, k, 1, 1);
        const auto state = c::makeStateSpace(A, B, s.readouts.state, {}, k, 1, d);
        // Predicted echo-only KLQE: measurement noise = readout residual + sensor noise.
        const double echoVariance = s.readouts.echoStd * s.readouts.echoStd;
        const double V = std::max((s.readouts.echoError * s.readouts.echoError
            + options.sensorNoise * options.sensorNoise) * echoVariance * dt, 1e-300);
        const auto lqe = c::designLQE(c::makeStateSpace(A, {}, s.readouts.echo, {}, k, 0, 1), W, {V});
        double total = 0;
        for (std::size_t i = 0; i < d; ++i) {
            double v = 0;
            for (std::size_t a = 0; a < k; ++a) for (std::size_t b = 0; b < k; ++b)
                v += s.readouts.state[i * k + a] * lqe.covariance[a * k + b] * s.readouts.state[i * k + b];
            const double sd = s.readouts.stateStd[i];
            const double e = std::sqrt(s.readouts.stateError[i] * s.readouts.stateError[i] + (sd > 0 ? v / (sd * sd) : 0.0));
            s.predictedStateErrors.push_back(e); total += e * e;
        }
        s.predictedStateError = std::sqrt(total / double(std::max<std::size_t>(d, 1)));
        s.observabilityBalance = minEigenvalueRatio(c::observabilityGramian(echo), k);
        s.pbhMargin = c::analyzeStructure(echo).pbhObservabilityMargin;
        // Output controllability of the state: D^-1 C_x W_c C_x^T D^-1.
        const auto Wc = c::controllabilityGramian(echo);
        std::vector<double> G(d * d, 0.0);
        for (std::size_t i = 0; i < d; ++i) for (std::size_t j = 0; j < d; ++j) {
            double v = 0;
            for (std::size_t a = 0; a < k; ++a) for (std::size_t b = 0; b < k; ++b)
                v += s.readouts.state[i * k + a] * Wc[a * k + b] * s.readouts.state[j * k + b];
            const double si = s.readouts.stateStd[i], sj = s.readouts.stateStd[j];
            G[i * d + j] = si > 0 && sj > 0 ? v / (si * sj) : 0.0;
        }
        double smallest = 0; (void)minEigenvalueRatio(G, d, &smallest);
        s.stateReach = std::sqrt(std::max(smallest, 0.0));
        const auto he = c::hankelSingularValues(echo), hs = c::hankelSingularValues(state);
        s.hsvEchoMax = he.front(); s.hsvEchoMin = he.back(); s.hsvStateMax = hs.front(); s.hsvStateMin = hs.back();
    } catch (const std::exception& error) { s.valid = false; s.failure = error.what(); }
    return s;
}

}  // namespace dynsyslab::workflow::postprocess_detail::mode_selection
