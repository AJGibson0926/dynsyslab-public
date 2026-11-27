// viewer_figure curves with their own source archive and clock (ViewerFigureJob.hpp).
//
// A curve may name "source" (another archive) and "clock" (the array its series are sampled on,
// default "time"). The figure's clock is "figure.clock": {"source": PATH (default the job source),
// "array": "time"}. Every curve is resampled onto that clock (core::TimeSeries::resample; linear,
// NaN outside the curve's own span; "interpolation" per curve picks linear, cubic, pchip or nearest)
// into one scratch archive, which the Viewer CLI then draws exactly like a single-source figure.
//
// Model-scale pages: "figure.fit_y_to": [curve index, ...] fits the ordinate range to those curves
// alone (within start/stop, 5% margin), so a small model response is not flattened by a large
// reference drawn on the same page.
#include "ViewerFigureJob.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/TimeSeriesOperations.hpp>

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace dynsyslab::workflow::detail {
namespace {

core::NpyArray onClock(const core::NpzRecord& record, const std::string& key, const std::string& clockKey,
                       const std::vector<double>& reference, const std::string& interpolation) {
    const auto& array = record.at(key);
    const auto clock = record.at(clockKey).doubles();
    const std::size_t width = array.shape.size() == 2 ? array.shape[1] : 1;
    if (array.shape.empty() || array.shape[0] != clock.size())
        throw std::invalid_argument("viewer_figure: series '" + key + "' is not sampled on '" + clockKey + "'");
    if (clock == reference) return array;  // already on the figure clock
    core::ResampleOptions options{core::interpolationFromName(interpolation), core::Extrapolation::nan};
    const auto resampled = core::TimeSeries(clock, array.doubles(), width).resample(reference, options);
    std::vector<std::size_t> shape{reference.size()};
    if (array.shape.size() == 2) shape.push_back(width);
    return core::NpyArray::fromDoubles(shape, resampled.valueVector());
}

}  // namespace

bool hasCurveSources(const json::Value& job) {
    if (!job.at("figure").at("clock").isNull()) return true;
    for (const auto& curve : job.at("figure").at("curves").array())
        if (!curve.at("source").isNull() || !curve.at("clock").isNull()) return true;
    return false;
}

json::Value assembleCurveSources(const json::Value& job, const std::filesystem::path& scratch) {
    auto resolved = job;
    const auto& clockSpec = job.at("figure").at("clock");
    const auto clockSource = clockSpec.at("source").string(job.at("source").string());
    const auto reference = core::readNpz(std::filesystem::u8path(clockSource))
                               .at(clockSpec.at("array").string("time")).doubles();
    if (reference.size() < 2) throw std::invalid_argument("viewer_figure: the figure clock needs two samples");
    core::NpzRecord out;
    out.arrays["time"] = core::NpyArray::fromDoubles({reference.size()}, reference);
    auto& curves = resolved.at("figure").at("curves").array();
    for (std::size_t k = 0; k < curves.size(); ++k) {
        auto& curve = curves[k];
        const auto source = core::readNpz(std::filesystem::u8path(curve.at("source").string(job.at("source").string())));
        const auto clockKey = curve.at("clock").string("time");
        const auto interpolation = curve.at("interpolation").string("linear");
        std::stringstream keys(curve.at("series").string());
        std::string key, renamed;
        while (std::getline(keys, key, ',')) {
            const auto name = "curve" + std::to_string(k) + "_" + key;
            out.arrays[name] = onClock(source, key, clockKey, reference, interpolation);
            renamed += (renamed.empty() ? "" : ",") + name;
        }
        curve.at("series") = json::Value(renamed);
    }
    out.arrays["metadata_json"] = core::NpyArray::fromString(
        R"({"schema":"dynsyslab.derived/v1","method":"viewer-figure-curves","scratch":true})");
    core::writeNpzAtomic(scratch, out);
    resolved.at("source") = json::Value(scratch.u8string());
    return resolved;
}

json::Value fitOrdinate(const json::Value& job) {
    const auto& chosen = job.at("figure").at("fit_y_to").array();
    if (chosen.empty() || !job.at("figure").at("y_range").isNull()) return job;
    const auto source = core::readNpz(std::filesystem::u8path(job.at("source").string()));
    const auto time = source.at(job.at("figure").at("time_key").string("time")).doubles();
    const double start = job.at("figure").at("start").number(-HUGE_VAL), stop = job.at("figure").at("stop").number(HUGE_VAL);
    double low = HUGE_VAL, high = -HUGE_VAL;
    const auto& curves = job.at("figure").at("curves").array();
    for (const auto& index : chosen) {
        const auto k = static_cast<std::size_t>(index.number());
        if (k >= curves.size()) throw std::invalid_argument("viewer_figure: fit_y_to names a curve that does not exist");
        auto key = curves[k].at("series").string();
        if (const auto comma = key.find(','); comma != std::string::npos) key = key.substr(comma + 1);
        const auto& array = source.at(key);
        const std::size_t width = array.shape.size() == 2 ? array.shape[1] : 1;
        const auto values = array.doubles();
        for (std::size_t r = 0; r < time.size() && r * width < values.size(); ++r) {
            if (time[r] < start || time[r] > stop) continue;
            for (std::size_t c = 0; c < width; ++c)
                if (std::isfinite(values[r * width + c])) { low = std::min(low, values[r * width + c]); high = std::max(high, values[r * width + c]); }
        }
    }
    if (!(high >= low)) throw std::invalid_argument("viewer_figure: the fit_y_to curves have no finite samples in the window");
    const double margin = 0.05 * std::max(high - low, 1e-12 * std::max(1.0, std::abs(high)));
    auto resolved = job;
    resolved.at("figure").at("y_range") = json::Value(json::Value::Array{json::Value(low - margin), json::Value(high + margin)});
    return resolved;
}

}  // namespace dynsyslab::workflow::detail
