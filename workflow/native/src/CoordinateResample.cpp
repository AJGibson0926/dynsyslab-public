#include "CoordinateChangeKinds.hpp"

#include <dynsyslab/core/TimeSeriesOperations.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

// coordinate_transform kind "resample": the source (a physical record or earlier coordinates) on
// a new clock, by core::resample with linear, cubic (Catmull-Rom Hermite), monotone PCHIP
// (Fritsch & Carlson 1980; Fritsch & Butland 1984) or nearest interpolation. It replaces the
// Python in-memory postprocessing/resampling.py, so non-uniform records can feed fixed-lag
// methods; the published record names its source and the method.

namespace dynsyslab::workflow::postprocess_detail {

dynsyslab::core::CoordinateSeries resampledCoordinates(const SourceData& source, const json::Value& settings,
                                                       json::Value& meta) {
    namespace c = dynsyslab::core;
    const auto n = source.samples();
    if (n < 2) throw std::invalid_argument("resample needs at least two source samples");
    const auto& time = source.time();
    std::vector<double> target;
    for (const auto& value : settings.at("times").array()) target.push_back(value.number());
    const double step = settings.at("time_step").number(0.0);
    if (!target.empty() && step > 0.0) throw std::invalid_argument("resample takes times or time_step, not both");
    if (target.empty()) {
        if (!(step > 0.0)) throw std::invalid_argument("resample needs times or a positive time_step");
        const double start = settings.at("start").number(time.front()), stop = settings.at("stop").number(time.back());
        const auto count = static_cast<std::size_t>(std::floor((stop - start) / step + 1.0e-9));
        for (std::size_t k = 0; k <= count; ++k) target.push_back(start + static_cast<double>(k) * step);
    }
    if (target.size() < 2) throw std::invalid_argument("resample needs at least two target times");
    for (std::size_t k = 1; k < target.size(); ++k)
        if (!(target[k] > target[k - 1])) throw std::invalid_argument("resample times must be strictly increasing");
    const double slack = 1.0e-12 * std::max({1.0, std::abs(time.front()), std::abs(time.back())});
    if (target.front() < time.front() - slack || target.back() > time.back() + slack)
        throw std::invalid_argument("resample times must lie within the source time span");
    for (auto& value : target) value = std::clamp(value, time.front(), time.back());
    const auto method = c::interpolationFromName(settings.at("method").string("pchip"));
    const auto view = c::TimeSeries::view(time, source.values().data(), n, source.dimension());
    auto series = view.resample(target, {method, c::Extrapolation::error});
    meta.object()["resampling"] = json::Value::Object{
        {"method", std::string(c::interpolationName(method))},
        {"raw_samples", static_cast<double>(n)},
        {"output_samples", static_cast<double>(target.size())},
        {"time_step", step},
        {"source_uniform_time_steps", view.uniformTimeSteps()}};
    meta.object()["uniform_time_steps"] = series.uniformTimeSteps();
    // Coordinate names and labels carry over: resampling keeps the components.
    const auto original = source.metadataJson.empty() ? json::Value{} : json::parse(source.metadataJson);
    for (const auto key : {"coordinate_names", "coordinate_labels", "coordinate_definition", "coordinate_provider",
                           "value_space", "coordinate_transform_identity"})
        if (!original.at(key).isNull() && settings.at("names").isNull()) meta.object()[key] = original.at(key);
    // A resampled physical record stays a physical record (state coordinates of its system).
    if (original.at("coordinate_provider").isNull()) {
        meta.object().erase("coordinate_provider");
        meta.object().erase("value_space");
    }
    return {std::move(series), std::nullopt};
}

}  // namespace dynsyslab::workflow::postprocess_detail
