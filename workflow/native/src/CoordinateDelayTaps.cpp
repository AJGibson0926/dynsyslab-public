#include "CoordinateChangeKinds.hpp"

#include <dynsyslab/core/DelaySampling.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

// coordinate_transform kind "delay_taps": delay coordinates at physical-time lags (core/DelaySampling.hpp).
//
//   {"kind": "delay_taps", "lags": [0, 0.15, 0.3],     physical lags (or "lag_samples": fractional ticks)
//    "columns": [0, 1],                                 source columns (default all)
//    "jet": [0, 1, 2]}                                  or: value, rate, acceleration columns, each tap
//                                                       read as a C^2 quintic delay jet (y, y', y'')
//
// The coordinates at sample t_k are y_c(t_k - tau_l), lag-major (all channels of the first lag, then
// the next), for every sample whose whole causal window lies in the record. Taps between samples are
// interpolated (linear, or the quintic Hermite jet), so lags need not be multiples of the sampling
// step and the window can be stated in physical time.

namespace dynsyslab::workflow::postprocess_detail {

dynsyslab::core::CoordinateSeries delayTapCoordinates(const SourceData& source, const json::Value& settings, json::Value& meta) {
    namespace c = dynsyslab::core;
    const auto n = source.samples(), d = source.dimension();
    if (n < 2) throw std::invalid_argument("delay_taps needs at least two source samples");
    std::vector<double> lags;
    for (const auto& value : settings.at("lags").array()) lags.push_back(value.number());
    if (lags.empty()) {
        const double dt = (source.time().back() - source.time().front()) / double(n - 1);
        for (const auto& value : settings.at("lag_samples").array()) lags.push_back(value.number() * dt);
    }
    if (lags.empty()) throw std::invalid_argument("delay_taps needs lags (physical time) or lag_samples");
    const double longest = *std::max_element(lags.begin(), lags.end());
    std::vector<double> endpoints;
    for (std::size_t k = 0; k < n; ++k)
        if (source.time()[k] - longest >= source.time().front() - 1e-12 * std::max(1.0, std::abs(source.time()[k])))
            endpoints.push_back(source.time()[k]);
    if (endpoints.size() < 2) throw std::invalid_argument("delay_taps: the longest lag leaves fewer than two complete windows");
    const auto view = c::TimeSeries::view(source.time(), source.values().data(), n, d);
    std::vector<std::string> names;
    const auto channelName = [&](std::size_t column) {
        return column < source.record.channels.size() && !source.record.channels[column].name.empty()
            ? source.record.channels[column].name : "q" + std::to_string(column);
    };
    const auto tag = [&](double lag) {
        std::string text = std::to_string(lag);
        text.erase(text.find_last_not_of('0') + 1);
        if (!text.empty() && text.back() == '.') text.pop_back();
        return "(t-" + text + ")";
    };
    std::vector<double> values;
    std::size_t width = 0;
    if (!settings.at("jet").isNull()) {
        std::vector<std::size_t> columns;
        for (const auto& value : settings.at("jet").array()) columns.push_back(std::size_t(value.number()));
        if (columns.size() != 3) throw std::invalid_argument("delay_taps jet needs value, rate and acceleration columns");
        std::vector<double> jet(n * 3);
        for (std::size_t k = 0; k < n; ++k)
            for (std::size_t i = 0; i < 3; ++i) {
                if (columns[i] >= d) throw std::invalid_argument("delay_taps jet column out of range");
                jet[k * 3 + i] = source.values()[k * d + columns[i]];
            }
        values = c::sampleDelayJet(c::TimeSeries(source.time(), std::move(jet), 3), endpoints, lags);
        width = 3 * lags.size();
        for (const double lag : lags)
            for (const char* part : {"", "'", "''"}) names.push_back(channelName(columns[0]) + part + tag(lag));
    } else {
        std::vector<std::size_t> columns;
        for (const auto& value : settings.at("columns").array()) columns.push_back(std::size_t(value.number()));
        if (columns.empty()) for (std::size_t c0 = 0; c0 < d; ++c0) columns.push_back(c0);
        const auto all = c::delaySamples(view, endpoints, lags);  // endpoints x lags x d
        for (std::size_t e = 0; e < endpoints.size(); ++e)
            for (std::size_t l = 0; l < lags.size(); ++l)
                for (const auto column : columns) {
                    if (column >= d) throw std::invalid_argument("delay_taps column out of range");
                    values.push_back(all[(e * lags.size() + l) * d + column]);
                }
        width = columns.size() * lags.size();
        for (const double lag : lags) for (const auto column : columns) names.push_back(channelName(column) + tag(lag));
    }
    json::Value::Array lagList;
    for (const double lag : lags) lagList.emplace_back(lag);
    meta.object()["delay_taps"] = json::Value::Object{{"lags", lagList}, {"interpolation", settings.at("jet").isNull() ? "linear" : "quintic-hermite-jet"},
        {"windows", double(endpoints.size())}, {"causal", true}};
    if (settings.at("names").isNull()) meta.object()["coordinate_names"] = json::parse(jsonStringArray(names));
    meta.object()["value_space"] = "history-coordinates";
    return {c::TimeSeries(std::move(endpoints), std::move(values), width), std::nullopt};
}

}  // namespace dynsyslab::workflow::postprocess_detail
