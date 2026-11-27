// Field-model definitions, shared regression and rollouts (FieldModels.hpp).
#include "FieldModelsInternal.hpp"

#include <dynsyslab/core/ModelLibraries.hpp>
#include <dynsyslab/core/SindyDynamics.hpp>
#include <dynsyslab/core/TimeSeries.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <optional>
#include <stdexcept>

namespace dynsyslab::core {
namespace detail {

json::Value numbersJson(const std::vector<double>& v) {
    json::Value::Array out;
    for (const auto x : v) out.emplace_back(x);
    return json::Value(std::move(out));
}

std::vector<double> jsonNumbers(const json::Value& v) {
    std::vector<double> out;
    for (const auto& x : v.array()) out.push_back(x.number());
    return out;
}

std::vector<std::pair<double, std::size_t>> nearest(const std::vector<double>& points, std::size_t d, const double* q, std::size_t k) {
    const auto rows = points.size() / d;
    std::vector<std::pair<double, std::size_t>> all(rows);
    for (std::size_t r = 0; r < rows; ++r) {
        double s = 0;
        for (std::size_t j = 0; j < d; ++j) s += (points[r * d + j] - q[j]) * (points[r * d + j] - q[j]);
        all[r] = {std::sqrt(s), r};
    }
    k = std::min(k, rows);
    std::partial_sort(all.begin(), all.begin() + std::ptrdiff_t(k), all.end());
    all.resize(k);
    return all;
}

std::vector<double> wendlandWeights(const std::vector<std::pair<double, std::size_t>>& hits) {
    const auto k = hits.size() - 1;
    const double radius = std::max(hits.back().first, 1e-12);
    std::vector<double> w(k);
    bool outside = true;
    for (std::size_t i = 0; i < k; ++i) {
        const double r = std::min(hits[i].first / radius, 1.0);
        outside = outside && r >= 1 - 1e-6;
        w[i] = std::pow(1 - r, 4) * (1 + 4 * r);
    }
    if (outside) std::fill(w.begin(), w.end(), 1.0);
    double s = 0;
    for (const auto v : w) s += v;
    for (auto& v : w) v /= s;
    return w;
}

std::vector<double> choleskySolve(std::vector<double> A, std::size_t n, std::vector<double> b, std::size_t m) {
    for (std::size_t j = 0; j < n; ++j) {
        double s = A[j * n + j];
        for (std::size_t k = 0; k < j; ++k) s -= A[j * n + k] * A[j * n + k];
        if (!(s > 0)) throw std::runtime_error("field regression: the normal equations are not positive definite");
        A[j * n + j] = std::sqrt(s);
        for (std::size_t i = j + 1; i < n; ++i) {
            double t = A[i * n + j];
            for (std::size_t k = 0; k < j; ++k) t -= A[i * n + k] * A[j * n + k];
            A[i * n + j] = t / A[j * n + j];
        }
    }
    for (std::size_t c = 0; c < m; ++c) {
        for (std::size_t i = 0; i < n; ++i) {
            double t = b[i * m + c];
            for (std::size_t k = 0; k < i; ++k) t -= A[i * n + k] * b[k * m + c];
            b[i * m + c] = t / A[i * n + i];
        }
        for (std::size_t i = n; i-- > 0;) {
            double t = b[i * m + c];
            for (std::size_t k = i + 1; k < n; ++k) t -= A[k * n + i] * b[k * m + c];
            b[i * m + c] = t / A[i * n + i];
        }
    }
    return b;
}

std::vector<double> trapezoidWeights(const std::vector<double>& t) {
    std::vector<double> w(t.size(), 0);
    for (std::size_t k = 1; k < t.size(); ++k) { const double h = (t[k] - t[k - 1]) / 2; w[k - 1] += h; w[k] += h; }
    return w;
}

}  // namespace detail

std::vector<double> thresholdedLeastSquares(const std::vector<double>& gram, std::size_t nf, std::size_t outputs, double ridge,
                                            double threshold) {
    const auto N = nf + outputs;
    if (gram.size() != N * N || !(ridge > 0) || threshold < 0) throw std::invalid_argument("thresholded least squares: moments, ridge or threshold");
    std::vector<double> norm(N);
    for (std::size_t i = 0; i < N; ++i) norm[i] = std::sqrt(std::max(gram[i * N + i], 1e-60));
    double top = 0;
    for (std::size_t i = 0; i < nf; ++i) top = std::max(top, norm[i]);
    std::vector<double> out(nf * outputs, 0);
    for (std::size_t j = 0; j < outputs; ++j) {
        std::vector<std::size_t> active;
        for (std::size_t i = 0; i < nf; ++i) if (norm[i] > 1e-10 * top) active.push_back(i);
        while (!active.empty()) {
            const auto a = active.size();
            std::vector<double> A(a * a), b(a);
            for (std::size_t r = 0; r < a; ++r) {
                for (std::size_t c = 0; c < a; ++c) A[r * a + c] = gram[active[r] * N + active[c]] / (norm[active[r]] * norm[active[c]]);
                A[r * a + r] += ridge;
                b[r] = gram[active[r] * N + nf + j] / (norm[active[r]] * norm[nf + j]);
            }
            const auto w = detail::choleskySolve(A, a, b, 1);
            std::vector<std::size_t> keep;
            for (std::size_t r = 0; r < a; ++r) if (std::abs(w[r]) >= threshold) keep.push_back(active[r]);
            if (keep.size() == a) {
                for (std::size_t r = 0; r < a; ++r) out[active[r] * outputs + j] = w[r] * norm[nf + j] / norm[active[r]];
                break;
            }
            active = std::move(keep);
        }
    }
    return out;
}

namespace {

Observable named(std::size_t inputs, std::size_t outputs, Observable::Evaluator f, const std::string& kind) {
    return Observable({kind, kind + " field model", "", "model", {}, false}, inputs, outputs, std::move(f));
}

}  // namespace

Observable fieldFromDefinition(const json::Value& d) {
    const auto kind = d.at("kind").string();
    if (kind == "polynomial_model") {
        const auto model = std::make_shared<SindyDynamicsModel>(SindyDynamicsModel::fromDefinition(d.at("model")));
        return named(model->dimension(), model->outputs(), [model](const double* x) {
            const auto v = model->evaluate(x);
            return Observable::Values(v.begin(), v.end());
        }, kind);
    }
    if (kind == "sum") {
        std::vector<Observable> parts;
        for (const auto& f : d.at("fields").array()) parts.push_back(fieldFromDefinition(f));
        if (parts.empty()) throw std::invalid_argument("sum field without fields");
        for (const auto& p : parts)
            if (p.stateDimension() != parts[0].stateDimension() || p.outputDimension() != parts[0].outputDimension())
                throw std::invalid_argument("sum field: parts differ in dimension");
        return named(parts[0].stateDimension(), parts[0].outputDimension(), [parts](const double* x) {
            auto v = parts[0].evaluate(x);
            for (std::size_t k = 1; k < parts.size(); ++k) { const auto w = parts[k].evaluate(x); for (std::size_t i = 0; i < v.size(); ++i) v[i] += w[i]; }
            return v;
        }, kind);
    }
    if (kind == "compact_radial") return detail::compactRadialField(d);
    if (kind == "local_tangent") return detail::localTangentField(d);
    if (kind == "paired_phase") return detail::pairedPhaseField(d);
    if (kind == "second_order") {
        const auto acceleration = fieldFromDefinition(d.at("acceleration"));
        const auto n = acceleration.stateDimension();
        if (n % 2 || acceleration.outputDimension() != n / 2) throw std::invalid_argument("second_order: acceleration of value/rate pairs");
        return named(n, n, [acceleration, n](const double* x) {
            const auto a = acceleration.evaluateReal(x);
            Observable::Values out(n);
            for (std::size_t i = 0; i < n / 2; ++i) { out[2 * i] = x[2 * i + 1]; out[2 * i + 1] = a[i]; }
            return out;
        }, kind);
    }
    if (kind == "pair_network") {
        const auto lags = detail::jsonNumbers(d.at("lags")), weights = detail::jsonNumbers(d.at("weights"));
        const auto n = lags.size();
        std::vector<Observable> fields;
        for (const auto& f : d.at("fields").array()) fields.push_back(fieldFromDefinition(f));
        std::vector<std::array<std::size_t, 3>> pairs;
        for (const auto& p : d.at("pairs").array())
            pairs.push_back({std::size_t(p.array().at(0).number()), std::size_t(p.array().at(1).number()), std::size_t(p.array().at(2).number())});
        return named(2 * n, 2 * n, [n, weights, fields, pairs](const double* a) {
            Observable::Values out(2 * n, 0.0);
            for (std::size_t i = 0; i < n; ++i) out[2 * i] = a[2 * i + 1];
            for (const auto& [i, j, k] : pairs) {
                const double z[4] = {a[2 * i], a[2 * i + 1], a[2 * j], a[2 * j + 1]};
                const auto v = fields.at(k).evaluateReal(z);
                out[2 * i + 1] += weights[i * n + j] * v.at(0);
                out[2 * j + 1] += weights[j * n + i] * v.at(1);
            }
            return out;
        }, kind);
    }
    if (kind == "interior") {
        const auto outer = fieldFromDefinition(d.at("outer")), drive = fieldFromDefinition(d.at("drive"));
        const bool feedback = !d.at("feedback").isNull();
        const auto correction = feedback ? fieldFromDefinition(d.at("feedback")) : outer;
        const auto gamma = detail::jsonNumbers(d.at("gamma"));
        const auto inner = drive.outputDimension(), n = 2 * (inner + 2);
        return named(n, n, [outer, drive, correction, feedback, gamma, inner, n](const double* z) {
            const double ends[4] = {z[0], z[1], z[n - 2], z[n - 1]};
            auto v = outer.evaluateReal(ends);
            if (feedback) { const auto c = correction.evaluateReal(z); v[0] += c[0]; v[1] += c[1]; }
            const auto f = drive.evaluateReal(ends);
            Observable::Values out(n);
            for (std::size_t i = 0; i < n / 2; ++i) out[2 * i] = z[2 * i + 1];
            out[1] = v[0]; out[n - 1] = v[1];
            for (std::size_t k = 0; k < inner; ++k) {
                const double g = gamma.size() == 1 ? gamma[0] : gamma.at(k);
                out[3 + 2 * k] = f[k] - 2 * g * z[3 + 2 * k] - g * g * z[2 + 2 * k];
            }
            return out;
        }, kind);
    }
    throw std::invalid_argument("unknown field model kind " + kind);
}

Rollout rolloutField(const Observable& f, const std::vector<double>& x0, const std::vector<double>& t, double maximumStep) {
    const auto n = x0.size();
    if (f.stateDimension() != n || f.outputDimension() != n || t.empty()) throw std::invalid_argument("field rollout: a vector field, a state and a clock");
    Rollout out;
    auto x = x0;
    const auto rate = [&](const std::vector<double>& s) { return f.evaluateReal(s.data()); };
    for (std::size_t k = 0; k < t.size(); ++k) {
        if (k) {
            const double span = t[k] - t[k - 1];
            const auto steps = maximumStep > 0 ? std::max<std::size_t>(1, std::size_t(std::ceil(std::abs(span) / maximumStep))) : 4;
            const double h = span / double(steps);
            for (std::size_t s = 0; s < steps; ++s) {
                const auto k1 = rate(x);
                std::vector<double> y(n);
                for (std::size_t i = 0; i < n; ++i) y[i] = x[i] + h / 2 * k1[i];
                const auto k2 = rate(y);
                for (std::size_t i = 0; i < n; ++i) y[i] = x[i] + h / 2 * k2[i];
                const auto k3 = rate(y);
                for (std::size_t i = 0; i < n; ++i) y[i] = x[i] + h * k3[i];
                const auto k4 = rate(y);
                for (std::size_t i = 0; i < n; ++i) x[i] += h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
            }
            if (!std::all_of(x.begin(), x.end(), [](double v) { return std::isfinite(v); })) return out;
        }
        out.time.push_back(t[k]);
        out.state.insert(out.state.end(), x.begin(), x.end());
    }
    out.complete = true;
    return out;
}

json::Value fitMap(const std::vector<FieldRecord>& records, const MapFitOptions& o, const json::Value& baseline) {
    if (records.empty()) throw std::invalid_argument("map fit: no records");
    const auto d = records[0].dimension, m = records[0].outputs;
    std::vector<double> all, time;
    for (const auto& r : records) {
        if (r.dimension != d || r.outputs != m || r.state.size() != r.samples() * d || r.label.size() != r.samples() * m)
            throw std::invalid_argument("map fit: record dimensions disagree");
        all.insert(all.end(), r.state.begin(), r.state.end());
        for (std::size_t k = 0; k < r.samples(); ++k) time.push_back(double(time.size()));
    }
    ModelLibraryOptions lo;
    lo.degree = o.degree; lo.centers = o.centers; lo.harmonics = o.harmonics; lo.rationalWidths = o.rationalWidths;
    const auto library = modelLibrary(o.family, TimeSeries(time, all, d), lo);
    const auto F = library.outputDimension(), N = F + m;
    const auto base = baseline.isNull() ? std::optional<Observable>{} : std::optional<Observable>{fieldFromDefinition(baseline)};
    std::vector<double> gram(N * N, 0);
    for (const auto& r : records) {
        const auto w = o.timeWeighted ? detail::trapezoidWeights(r.time) : std::vector<double>(r.samples(), 1.0);
        for (std::size_t k = 0; k < r.samples(); ++k) {
            auto row = library.evaluateReal(r.state.data() + k * d);
            for (std::size_t j = 0; j < m; ++j) row.push_back(r.label[k * m + j]);
            if (base) { const auto b = base->evaluateReal(r.state.data() + k * d); for (std::size_t j = 0; j < m; ++j) row[F + j] -= b.at(j); }
            const double s = o.timeWeighted ? std::sqrt(w[k]) : 1.0;
            for (auto& v : row) v *= s;
            for (std::size_t a = 0; a < N; ++a) for (std::size_t b = 0; b < N; ++b) gram[a * N + b] += row[a] * row[b];
        }
    }
    const SindyDynamicsModel model(library, thresholdedLeastSquares(gram, F, m, o.ridge, o.threshold), m);
    json::Value fitted = json::Value::Object{{"kind", "polynomial_model"}, {"model", model.definition()}, {"family", o.family}};
    if (!base) return fitted;
    return json::Value::Object{{"kind", "sum"}, {"fields", json::Value::Array{baseline, fitted}}};
}

}  // namespace dynsyslab::core
