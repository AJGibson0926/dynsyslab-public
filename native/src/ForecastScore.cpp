// Forecast scoring on independent clocks, prominent-peak matching and the kinematic chain-rule check
// (ForecastScore.hpp).
#include <dynsyslab/core/ForecastScore.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::core {

std::vector<std::size_t> prominentPeaks(const std::vector<double>& x, double prominence) {
    if (!(prominence >= 0) || !std::isfinite(prominence)) throw std::invalid_argument("peak prominence must be finite and nonnegative");
    std::vector<std::size_t> peaks;
    const std::size_t n = x.size();
    for (std::size_t i = 1; i + 1 < n; ++i) {
        if (!(x[i] > x[i-1])) continue;
        std::size_t j = i;  // a plateau's peak is its middle sample
        while (j + 1 < n && x[j+1] == x[i]) ++j;
        if (j + 1 >= n || !(x[j+1] < x[i])) continue;
        const std::size_t peak = (i + j) / 2;
        // Prominence: the peak's height above the higher of its two bases (the lowest point before a
        // higher sample, or the record's end, on each side).
        double left = x[i];
        for (std::size_t k = i; k-- > 0;) { if (x[k] > x[i]) break; left = std::min(left, x[k]); }
        double right = x[i];
        for (std::size_t k = j + 1; k < n; ++k) { if (x[k] > x[i]) break; right = std::min(right, x[k]); }
        if (x[i] - std::max(left, right) >= prominence) peaks.push_back(peak);
        i = j;
    }
    return peaks;
}

json::Value forecastScore(const TimeSeries& reference, const TimeSeries& prediction, const ForecastScoreOptions& options) {
    const std::size_t n = reference.sampleCount(), d = reference.componentCount(), m = prediction.sampleCount();
    if (n < 2 || prediction.componentCount() != d || options.channel >= d)
        throw std::invalid_argument("forecast and reference must share their channels");
    const double* t = reference.timeData();
    const double* pt = prediction.timeData();
    json::Value::Object out{{"requested", json::Value::Array{t[0], t[n-1]}}, {"reference_count", double(n)}, {"prediction_count", double(m)}};
    std::vector<std::size_t> shared;
    for (std::size_t k = 0; m && k < n; ++k) if (t[k] >= pt[0] && t[k] <= pt[m-1]) shared.push_back(k);
    out["sample_coverage"] = double(shared.size()) / double(n);
    out["complete"] = shared.size() == n;
    out["compared_count"] = double(shared.size());
    if (shared.empty()) { out["coverage"] = 0.0; return out; }
    const double covered = std::max(0.0, std::min(pt[m-1], t[n-1]) - std::max(pt[0], t[0]));
    out["coverage"] = covered / (t[n-1] - t[0]);
    out["compared"] = json::Value::Array{t[shared.front()], t[shared.back()]};
    // Errors on the shared reference samples, the prediction linear between its own samples.
    std::vector<double> weight(shared.size(), 1.0);
    if (options.timeWeighted && shared.size() > 1) {
        std::fill(weight.begin(), weight.end(), 0.0);
        for (std::size_t s = 1; s < shared.size(); ++s) {
            const double h = t[shared[s]] - t[shared[s-1]];
            weight[s-1] += h / 2; weight[s] += h / 2;
        }
    }
    double total = 0;
    for (const double w : weight) total += w;
    json::Value::Array rms, maximum, relative;
    for (std::size_t c = 0; c < d; ++c) {
        double squared = 0, worst = 0, energy = 0;
        for (std::size_t s = 0; s < shared.size(); ++s) {
            const auto k = shared[s];
            const auto at = std::size_t(std::lower_bound(pt, pt + m, t[k]) - pt);
            double value = prediction.value(std::min(at, m - 1), c);
            if (at > 0 && at < m && pt[at] != t[k]) {
                const double a = (t[k] - pt[at-1]) / (pt[at] - pt[at-1]);
                value = prediction.value(at - 1, c) + a * (prediction.value(at, c) - prediction.value(at - 1, c));
            }
            const double error = value - reference.value(k, c);
            squared += weight[s] * error * error;
            energy += weight[s] * reference.value(k, c) * reference.value(k, c);
            worst = std::max(worst, std::abs(error));
        }
        rms.emplace_back(std::sqrt(squared / total));
        maximum.emplace_back(worst);
        relative.emplace_back(energy > 0 ? json::Value(std::sqrt(squared / energy)) : json::Value());
    }
    out["errors"] = json::Value::Object{{"rmse", rms}, {"maximum_absolute", maximum}, {"relative_rms", relative},
                                        {"time_weighted", options.timeWeighted}};
    if (options.prominence) {
        std::vector<double> x(n), p(m);
        for (std::size_t k = 0; k < n; ++k) x[k] = reference.value(k, options.channel);
        for (std::size_t k = 0; k < m; ++k) p[k] = prediction.value(k, options.channel);
        const auto truth = prominentPeaks(x, *options.prominence), peaks = prominentPeaks(p, *options.prominence);
        json::Value::Array matches;
        for (std::size_t j = 0; j < truth.size(); ++j) {
            const double left = j ? (t[truth[j-1]] + t[truth[j]]) / 2 : t[0];
            const double right = j + 1 < truth.size() ? (t[truth[j]] + t[truth[j+1]]) / 2 : t[n-1];
            json::Value::Object item{{"reference_time", t[truth[j]]}, {"reference_value", x[truth[j]]},
                                     {"full_interval", m && pt[0] <= left && pt[m-1] >= right}};
            std::size_t count = 0, best = m;
            for (const auto k : peaks)
                if (pt[k] >= left && pt[k] < right) { ++count; if (best == m || p[k] > p[best]) best = k; }
            item["count"] = double(count);
            if (best < m) {
                item["time"] = pt[best];
                item["value"] = p[best];
                item["time_error"] = pt[best] - t[truth[j]];
                item["relative_height"] = x[truth[j]] != 0 ? json::Value(p[best] / x[truth[j]] - 1) : json::Value();
            }
            matches.emplace_back(json::Value(std::move(item)));
        }
        std::size_t inside = 0;
        for (const auto k : peaks) inside += pt[k] >= t[0] && pt[k] <= t[n-1];
        out["peaks"] = json::Value::Object{{"reference_count", double(truth.size())}, {"forecast_count", double(inside)},
                                           {"prominence", *options.prominence}, {"matches", matches}};
    }
    return out;
}

double kinematicMismatch(const SindyDynamicsModel& model, const TimeSeries& q, std::size_t value, std::size_t rate) {
    if (!model.isVectorField() || value >= model.dimension() || rate >= model.dimension() || q.componentCount() != model.dimension())
        throw std::invalid_argument("the kinematic check needs a vector field with value and rate coordinates");
    double gap = 0, scale = 0;
    for (std::size_t k = 0; k < q.sampleCount(); ++k) {
        const auto f = model.evaluate(q.valuesData() + k * q.componentCount());
        gap += std::pow(f[value] - q.value(k, rate), 2);
        scale += std::pow(q.value(k, rate), 2);
    }
    return scale > 0 ? std::sqrt(gap / scale) : 0.0;
}

}  // namespace dynsyslab::core
