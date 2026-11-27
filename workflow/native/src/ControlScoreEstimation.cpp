// control_score "estimation" scores (ControlScoreEstimation.hpp): state-estimator and shadow-controller
// errors on the estimator clock. The truth is the solved state interpolated linearly onto that clock.
// Innovation whiteness is the normalized sample autocorrelation at time lags (core::autocorrelation;
// Box & Jenkins, Time Series Analysis, 1970; Mehra, IEEE Trans. Autom. Control 15(2), 175, 1970,
// for innovation-based filter tuning).
#include "ControlScoreEstimation.hpp"

#include <dynsyslab/core/SignalMetrics.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace {
namespace c = dynsyslab::core;
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

std::vector<std::size_t> components(const json::Value& estimation) {
    std::vector<std::size_t> out;
    for (const auto& item : estimation.at("components").array()) out.push_back(static_cast<std::size_t>(item.number()));
    return out;
}

std::string lagName(double lag) {
    char text[32];
    std::snprintf(text, sizeof(text), "innovation_autocorrelation_%g", lag);
    return text;
}

double interpolate(const std::vector<double>& xp, const std::vector<double>& fp, double x) {
    if (x <= xp.front()) return fp.front();
    if (x >= xp.back()) return fp.back();
    const auto j = static_cast<std::size_t>(std::upper_bound(xp.begin(), xp.end(), x) - xp.begin()) - 1;
    return fp[j] + (fp[j + 1] - fp[j]) * (x - xp[j]) / (xp[j + 1] - xp[j]);
}

double rms(const std::vector<double>& v) {
    if (v.empty()) return kNaN;
    double s = 0.0; for (const double x : v) s += x * x;
    return std::sqrt(s / double(v.size()));
}
double mean(const std::vector<double>& v) {
    if (v.empty()) return kNaN;
    double s = 0.0; for (const double x : v) s += x; return s / double(v.size());
}
double deviation(const std::vector<double>& v) {
    const double m = mean(v); std::vector<double> d; for (const double x : v) d.push_back(x - m); return rms(d);
}

}  // namespace

std::vector<std::string> estimationMetricNames(const json::Value& estimation) {
    std::vector<std::string> names;
    if (estimation.isNull()) return names;
    for (const auto k : components(estimation))
        for (const char* stem : {"estimate_error_rms_", "estimate_bias_", "state_std_", "estimate_relative_error_"})
            names.push_back(stem + std::to_string(k));
    names.push_back("innovation_relative");
    for (const auto& lag : estimation.at("innovation_lags").array()) names.push_back(lagName(lag.number()));
    if (estimation.at("shadow_force").boolean(false))
        for (const char* name : {"shadow_force_error_rms", "shadow_force_rms", "shadow_force_relative_error", "shadow_force_correlation"})
            names.emplace_back(name);
    return names;
}

std::map<std::string, double> scoreEstimation(const core::NpzRecord& record, double start, double stop,
                                              const json::Value& estimation) {
    std::map<std::string, double> out;
    if (estimation.isNull() || !record.contains("estimator_time")) return out;
    const auto clock = record.at("estimator_time").doubles();
    std::vector<std::size_t> rows;
    for (std::size_t k = 0; k < clock.size(); ++k) if (clock[k] >= start && clock[k] <= stop) rows.push_back(k);
    if (rows.size() < 2) return out;
    if (record.contains("state_estimate")) {
        const auto time = record.at("time").doubles();
        const auto values = record.at("values").doubles();
        const std::size_t dim = record.at("values").shape.size() > 1 ? record.at("values").shape[1] : 1;
        const auto estimate = record.at("state_estimate").doubles();
        const std::size_t edim = estimate.size() / clock.size();
        for (const auto k : components(estimation)) {
            if (k >= dim || k >= edim) throw std::invalid_argument("control_score estimation component is outside the state");
            std::vector<double> component(time.size()), truth, error;
            for (std::size_t r = 0; r < time.size(); ++r) component[r] = values[r * dim + k];
            for (const auto r : rows) {
                truth.push_back(interpolate(time, component, clock[r]));
                error.push_back(estimate[r * edim + k] - truth.back());
            }
            const auto id = std::to_string(k);
            out["estimate_error_rms_" + id] = rms(error);
            out["estimate_bias_" + id] = mean(error);
            out["state_std_" + id] = deviation(truth);
            out["estimate_relative_error_" + id] = rms(error) / deviation(truth);
        }
    }
    if (record.contains("innovation")) {
        const auto innovation = record.at("innovation").doubles();
        const std::size_t width = innovation.size() / clock.size();
        std::vector<double> window;
        for (const auto r : rows) window.push_back(innovation[r * width]);
        if (record.contains("measurement_values")) {
            const auto measured = record.at("measurement_values").doubles();
            const std::size_t mwidth = measured.size() / clock.size();
            std::vector<double> y; for (const auto r : rows) y.push_back(measured[r * mwidth]);
            out["innovation_relative"] = rms(window) / deviation(y);
        }
        const double dt = clock[rows[1]] - clock[rows[0]];
        for (const auto& lag : estimation.at("innovation_lags").array()) {
            const auto steps = std::max<std::size_t>(1, static_cast<std::size_t>(std::lround(lag.number() / dt)));
            out[lagName(lag.number())] = c::autocorrelation(window, steps);
        }
    }
    if (estimation.at("shadow_force").boolean(false) && record.contains("commanded_input") && record.contains("applied_input")) {
        const auto commanded = record.at("commanded_input").doubles(), applied = record.at("applied_input").doubles();
        std::vector<double> u, uk, difference;
        for (const auto r : rows) { u.push_back(applied.at(r)); uk.push_back(commanded.at(r)); difference.push_back(uk.back() - u.back()); }
        out["shadow_force_error_rms"] = rms(difference);
        out["shadow_force_rms"] = rms(u);
        out["shadow_force_relative_error"] = rms(difference) / rms(u);
        const double mu = mean(u), mk = mean(uk);
        double suk = 0, suu = 0, skk = 0;
        for (std::size_t i = 0; i < u.size(); ++i) { suk += (u[i] - mu) * (uk[i] - mk); suu += (u[i] - mu) * (u[i] - mu); skk += (uk[i] - mk) * (uk[i] - mk); }
        out["shadow_force_correlation"] = suu > 0 && skk > 0 ? suk / std::sqrt(suu * skk) : kNaN;
    }
    return out;
}

}  // namespace dynsyslab::workflow::postprocess_detail
