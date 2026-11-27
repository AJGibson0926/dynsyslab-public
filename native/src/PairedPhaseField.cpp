// Paired amplitude/phase fields (FieldModels.hpp). Ported from
// tools/experiments/miscellaneous_computes/paired.py. Each real pair (x_{2j}, x_{2j+1}) is a complex
// amplitude u_j; its rate is u_j' = (a_j + i w_j) u_j with radial and angular rates a_j, w_j linear in
// saturating amplitude features and, with `phase`, the relative phase of u_j to the harmonic
// u_0^{h_j} (phase locking at winding ratio h_j; Pikovsky, Rosenblum and Kurths 2001,
// "Synchronization", chapter 7). Fits use instantaneous rates and, optionally, weak-form integral
// rows over sample windows (Messenger and Bortz 2021, J. Comput. Phys. 443:110525).
#include "FieldModelsInternal.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <memory>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

constexpr double pi = 3.14159265358979323846;

std::size_t featureCount(bool phase) { return phase ? 12 : 6; }

void features(const double* b, std::size_t pairs, std::size_t j, int h, bool phase, double* f) {
    std::vector<std::complex<double>> u(pairs);
    std::vector<double> p(pairs);
    double E = 0;
    for (std::size_t k = 0; k < pairs; ++k) { u[k] = {b[2 * k], b[2 * k + 1]}; p[k] = std::norm(u[k]); E += p[k]; }
    E /= 2.0 * double(pairs);
    for (std::size_t k = 0; k < pairs; ++k) u[k] /= std::sqrt(p[k] + 0.008);
    const double q = p[j], p0 = p[0];
    f[0] = 1; f[1] = E / (E + 0.15); f[2] = E / (E + 1.5); f[3] = E / (E + 8); f[4] = p0 / (1 + p0); f[5] = q / (1 + q);
    if (!phase) return;
    const auto a = u[0];
    const auto base = std::pow(a, std::abs(h));
    const auto offset = u[j] * (h >= 0 ? std::conj(base) : base);
    const auto a2 = a * a;
    f[6] = a.real(); f[7] = a.imag(); f[8] = a2.real(); f[9] = a2.imag(); f[10] = offset.real(); f[11] = offset.imag();
}

struct Paired {
    std::vector<int> harmonics;
    bool phase{};
    std::vector<double> coef;  // pairs x features x 2
    [[nodiscard]] std::vector<double> operator()(const double* x) const {
        const auto pairs = harmonics.size(), F = featureCount(phase);
        std::vector<double> out(2 * pairs), f(F);
        for (std::size_t j = 0; j < pairs; ++j) {
            features(x, pairs, j, harmonics[j], phase, f.data());
            double a = 0, w = 0;
            for (std::size_t s = 0; s < F; ++s) { a += f[s] * coef[(j * F + s) * 2]; w += f[s] * coef[(j * F + s) * 2 + 1]; }
            const double u = x[2 * j], v = x[2 * j + 1];
            out[2 * j] = a * u - w * v; out[2 * j + 1] = w * u + a * v;
        }
        return out;
    }
};

std::vector<double> unwrapped(std::vector<double> angle) {
    for (std::size_t k = 1; k < angle.size(); ++k) {
        double d = angle[k] - angle[k - 1];
        d -= 2 * pi * std::round(d / (2 * pi));
        angle[k] = angle[k - 1] + d;
    }
    return angle;
}

}  // namespace

namespace detail {
Observable pairedPhaseField(const json::Value& d) {
    auto p = std::make_shared<Paired>();
    for (const auto& h : d.at("harmonics").array()) p->harmonics.push_back(int(h.number()));
    p->phase = d.at("phase").boolean();
    p->coef = jsonNumbers(d.at("coefficients"));
    const auto n = 2 * p->harmonics.size();
    if (!n || p->coef.size() != p->harmonics.size() * featureCount(p->phase) * 2)
        throw std::invalid_argument("paired_phase: two coefficient columns per pair");
    const std::shared_ptr<const Paired> field = p;
    return Observable({"paired_phase", "paired amplitude/phase field", "", "model", {}, false}, n, n,
                      [field](const double* x) { const auto y = (*field)(x); return Observable::Values(y.begin(), y.end()); });
}
}  // namespace detail

json::Value fitPairedPhaseField(const std::vector<FieldRecord>& records, const PairedPhaseOptions& o) {
    if (records.empty()) throw std::invalid_argument("paired_phase fit: no records");
    const auto n = records[0].dimension;
    if (n % 2 || !n || o.instantWeight < 0 || o.instantWeight > 1 || (o.instantWeight < 1 && o.weakWindows.empty()))
        throw std::invalid_argument("paired_phase fit: modal pairs and valid derivative/integral weights");
    for (const auto& r : records)
        if (r.dimension != n || r.outputs != n || r.samples() < 2) throw std::invalid_argument("paired_phase fit: records need rates of their own states");
    const auto pairs = n / 2;
    Paired model;
    model.phase = o.phase;
    model.harmonics = o.harmonics;
    if (model.harmonics.empty()) {
        const auto& r = records[0];
        std::vector<double> turns(pairs);
        for (std::size_t j = 0; j < pairs; ++j) {
            std::vector<double> angle(r.samples());
            for (std::size_t k = 0; k < r.samples(); ++k) angle[k] = std::atan2(r.state[k * n + 2 * j + 1], r.state[k * n + 2 * j]);
            const auto un = unwrapped(angle);
            turns[j] = un.back() - un.front();
        }
        if (std::abs(turns[0]) < 1e-8) throw std::invalid_argument("paired_phase fit: specify harmonics for a nonwinding reference pair");
        for (const auto t : turns) model.harmonics.push_back(std::clamp(int(std::lround(t / turns[0])), -24, 24));
    }
    if (model.harmonics.size() != pairs) throw std::invalid_argument("paired_phase fit: one winding number per pair");
    const auto F = featureCount(o.phase), N = F + 2;
    for (std::size_t j = 0; j < pairs; ++j) {
        std::vector<double> gram(N * N, 0), row(N);
        const auto add = [&](double weight) {
            for (std::size_t a = 0; a < N; ++a) for (std::size_t c = 0; c < N; ++c) gram[a * N + c] += weight * row[a] * row[c];
        };
        for (const auto& r : records) {
            const auto S = r.samples();
            std::vector<double> f(S * F), power(S), logPower(S), angle(S);
            for (std::size_t k = 0; k < S; ++k) {
                const double* b = r.state.data() + k * n;
                const double* v = r.label.data() + k * n;
                features(b, pairs, j, model.harmonics[j], o.phase, f.data() + k * F);
                const double x0 = b[2 * j], x1 = b[2 * j + 1];
                power[k] = x0 * x0 + x1 * x1;
                logPower[k] = 0.5 * std::log(power[k] + 1e-10);
                angle[k] = std::atan2(x1, x0);
                const double radius = std::sqrt(std::max(power[k], 1e-5));
                for (std::size_t s = 0; s < F; ++s) row[s] = f[k * F + s] * radius;
                row[F] = (x0 * v[2 * j] + x1 * v[2 * j + 1]) / radius;
                row[F + 1] = (x0 * v[2 * j + 1] - x1 * v[2 * j]) / radius;
                if (o.instantWeight > 0) add(o.instantWeight);
            }
            if (o.instantWeight >= 1) continue;
            angle = unwrapped(angle);
            std::vector<double> integral(S * F, 0);  // cumulative trapezoid of the features
            for (std::size_t k = 1; k < S; ++k)
                for (std::size_t s = 0; s < F; ++s)
                    integral[k * F + s] = integral[(k - 1) * F + s] + 0.5 * (r.time[k] - r.time[k - 1]) * (f[k * F + s] + f[(k - 1) * F + s]);
            for (const auto w : o.weakWindows) {
                if (w < 1 || w >= S) throw std::invalid_argument("paired_phase fit: invalid weak window length");
                for (std::size_t k = w; k < S; ++k) {
                    const double span = r.time[k] - r.time[k - w];
                    const double scale = std::pow(power[k] * power[k - w], 0.25) / span;
                    for (std::size_t s = 0; s < F; ++s) row[s] = (integral[k * F + s] - integral[(k - w) * F + s]) * scale;
                    row[F] = (logPower[k] - logPower[k - w]) * scale;
                    row[F + 1] = (angle[k] - angle[k - w]) * scale;
                    add((1 - o.instantWeight) / double(o.weakWindows.size()));
                }
            }
        }
        const auto c = thresholdedLeastSquares(gram, F, 2, o.ridge, o.threshold);
        model.coef.insert(model.coef.end(), c.begin(), c.end());
    }
    json::Value::Array harmonics;
    for (const auto h : model.harmonics) harmonics.emplace_back(double(h));
    return json::Value::Object{{"kind", "paired_phase"}, {"dimension", double(n)}, {"harmonics", harmonics}, {"phase", o.phase},
                               {"coefficients", detail::numbersJson(model.coef)}};
}

}  // namespace dynsyslab::core
