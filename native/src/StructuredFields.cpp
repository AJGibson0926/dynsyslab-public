// Structured second-order field fits (FieldModels.hpp). Ported from
// tools/experiments/miscellaneous_computes/structured.py: imposed kinematic rows with a fitted
// acceleration map (second-order SINDy; Brunton, Proctor and Kutz 2016), acceleration maps shared
// by lag separation across a tap network, and interior taps as critically damped responses to an
// endpoint drive, optionally with all-tap endpoint feedback.
#include "FieldModelsInternal.hpp"

#include <cmath>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

void requireRates(const std::vector<FieldRecord>& records, const char* what) {
    if (records.empty()) throw std::invalid_argument(std::string(what) + ": no records");
    const auto n = records[0].dimension;
    for (const auto& r : records)
        if (r.dimension != n || r.outputs != n || r.state.size() != r.samples() * n || r.label.size() != r.samples() * n)
            throw std::invalid_argument(std::string(what) + ": records need rates of their own states");
}

// A record whose state and label are chosen columns of another's (or computed per sample).
template <class State, class Label>
FieldRecord project(const FieldRecord& r, std::size_t d, std::size_t m, State state, Label label) {
    FieldRecord out{r.time, {}, {}, d, m};
    out.state.reserve(r.samples() * d); out.label.reserve(r.samples() * m);
    for (std::size_t k = 0; k < r.samples(); ++k) {
        const double* x = r.state.data() + k * r.dimension;
        const double* v = r.label.data() + k * r.outputs;
        for (std::size_t i = 0; i < d; ++i) out.state.push_back(state(x, i));
        for (std::size_t j = 0; j < m; ++j) out.label.push_back(label(x, v, j));
    }
    return out;
}

json::Value array(const std::vector<double>& v) { return detail::numbersJson(v); }

}  // namespace

json::Value fitSecondOrderField(const std::vector<FieldRecord>& records, const MapFitOptions& options) {
    requireRates(records, "second_order fit");
    const auto n = records[0].dimension;
    if (n % 2) throw std::invalid_argument("second_order fit: value/rate pairs");
    std::vector<FieldRecord> data;
    for (const auto& r : records)
        data.push_back(project(r, n, n / 2, [](const double* x, std::size_t i) { return x[i]; },
                               [](const double*, const double* v, std::size_t j) { return v[2 * j + 1]; }));
    return json::Value::Object{{"kind", "second_order"}, {"acceleration", fitMap(data, options)}};
}

json::Value fitPairNetwork(const std::vector<FieldRecord>& records, const std::vector<double>& lags, const std::vector<double>& weights,
                           const MapFitOptions& options, double tolerance) {
    requireRates(records, "pair_network fit");
    const auto n = lags.size();
    if (n < 2 || lags[0] != 0 || weights.size() != n * n || records[0].dimension != 2 * n || !std::isfinite(tolerance) || tolerance < 0)
        throw std::invalid_argument("pair_network fit: inconsistent design");
    for (std::size_t i = 0; i < n; ++i) {
        if (i && !(lags[i] > lags[i - 1])) throw std::invalid_argument("pair_network fit: lags increase from zero");
        double row = 0;
        for (std::size_t j = 0; j < n; ++j) {
            const double w = weights[i * n + j];
            if (!std::isfinite(w) || w < 0 || (i == j && w != 0)) throw std::invalid_argument("pair_network fit: weights nonnegative, zero diagonal");
            row += w;
        }
        if (std::abs(row - 1) > 1e-8) throw std::invalid_argument("pair_network fit: weight rows sum to one");
    }
    std::vector<double> gaps;
    std::vector<std::vector<FieldRecord>> groups;
    json::Value::Array pairs;
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = i + 1; j < n; ++j) {
            if (!(weights[i * n + j] || weights[j * n + i])) continue;
            const double gap = lags[j] - lags[i];
            std::size_t k = 0;
            while (k < gaps.size() && std::abs(gap - gaps[k]) > tolerance) ++k;
            if (k == gaps.size()) { gaps.push_back(gap); groups.emplace_back(); }
            pairs.push_back(json::Value::Array{double(i), double(j), double(k)});
            const std::size_t cols[4] = {2 * i, 2 * i + 1, 2 * j, 2 * j + 1}, rates[2] = {2 * i + 1, 2 * j + 1};
            for (const auto& r : records)
                groups[k].push_back(project(r, 4, 2, [&](const double* x, std::size_t c) { return x[cols[c]]; },
                                            [&](const double*, const double* v, std::size_t c) { return v[rates[c]]; }));
        }
    json::Value::Array fields;
    for (const auto& g : groups) fields.push_back(fitMap(g, options));
    return json::Value::Object{{"kind", "pair_network"}, {"lags", array(lags)}, {"weights", array(weights)},
                               {"pairs", pairs}, {"fields", fields}};
}

json::Value fitInteriorField(const std::vector<FieldRecord>& records, const InteriorOptions& o) {
    requireRates(records, "interior fit");
    const auto n = records[0].dimension, inner = n / 2 - 2;
    if (n < 6 || n % 2) throw std::invalid_argument("interior fit: two outer taps and at least one inner tap");
    if (o.gamma.empty() || (o.gamma.size() != 1 && o.gamma.size() != inner)) throw std::invalid_argument("interior fit: one damping rate or one per interior tap");
    for (const auto g : o.gamma) if (!std::isfinite(g) || g <= 0) throw std::invalid_argument("interior fit: damping must be positive");
    const auto gamma = [&](std::size_t k) { return o.gamma.size() == 1 ? o.gamma[0] : o.gamma[k]; };
    const auto end = [n](const double* x, std::size_t i) { const std::size_t c[4] = {0, 1, n - 2, n - 1}; return x[c[i]]; };
    std::vector<FieldRecord> outerData;
    for (const auto& r : records)
        outerData.push_back(project(r, 4, 2, end, [n](const double*, const double* v, std::size_t j) { return v[j ? n - 1 : 1]; }));
    const auto outer = fitMap(outerData, o.outer);
    json::Value feedback;
    if (o.feedback) {
        const auto outerField = fieldFromDefinition(outer);
        std::vector<FieldRecord> data;
        for (const auto& r : records)
            data.push_back(project(r, n, 2, [](const double* x, std::size_t i) { return x[i]; },
                                   [&](const double* x, const double* v, std::size_t j) {
                                       const double e[4] = {x[0], x[1], x[n - 2], x[n - 1]};
                                       return v[j ? n - 1 : 1] - outerField.evaluateReal(e)[j];
                                   }));
        feedback = fitMap(data, o.feedbackOptions);
    }
    json::Value drive;
    if (o.waveformDrive) {
        if (o.gamma.size() != 1) throw std::invalid_argument("interior fit: waveform fitting uses a common scalar damping rate");
        std::vector<DriveRecord> data;
        for (const auto& r : records) {
            DriveRecord d{r.time, {}, {}, {}, 4, inner};
            for (std::size_t k = 0; k < r.samples(); ++k) {
                const double* x = r.state.data() + k * n;
                for (std::size_t i = 0; i < 4; ++i) d.input.push_back(end(x, i));
                for (std::size_t i = 0; i < inner; ++i) { d.position.push_back(x[2 + 2 * i]); d.velocity.push_back(x[3 + 2 * i]); }
            }
            data.push_back(std::move(d));
        }
        drive = fitDampedDrive(data, o.gamma[0], o.drive);
    } else {
        std::vector<FieldRecord> data;
        for (const auto& r : records)
            data.push_back(project(r, 4, inner, end, [&](const double* x, const double* v, std::size_t k) {
                const double g = gamma(k);
                return v[3 + 2 * k] + 2 * g * x[3 + 2 * k] + g * g * x[2 + 2 * k];
            }));
        drive = fitMap(data, o.drive);
    }
    return json::Value::Object{{"kind", "interior"}, {"outer", outer}, {"drive", drive}, {"feedback", feedback}, {"gamma", array(o.gamma)}};
}

}  // namespace dynsyslab::core
