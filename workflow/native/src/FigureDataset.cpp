// Post-processing method "figure_dataset": a plotting table assembled from existing archives.
//
// Settings: {"inputs":{"run":PATH,...},
//            "arrays":[{"name":"radius","from":"run","array":"values","column":0,
//                       "onto":{"from":"base","array":"time"},"clock":"time","scale":1.0,
//                       "scale_before":1.0 (applied before interpolation; "scale" after),
//                       "offset":0.0 (added after "scale"),
//                       "keep_times":[t,...] (NaN except at the nearest samples, for point markers),
//                       "nan_after_end":false}, ...]}
// "columns":[0,1] keeps a matrix. "onto" interpolates linearly onto another array's clock, with
// the arithmetic of numpy.interp, and "nan_after_end" marks samples past a run's last time.
#include "PostprocessInternal.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Text.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace {

namespace c = dynsyslab::core;

std::vector<double> interpolate(const std::vector<double>& x, const std::vector<double>& xp,
                                const std::vector<double>& fp, bool nanAfterEnd) {
    if (xp.empty() || xp.size() != fp.size()) throw std::invalid_argument("figure_dataset: clock and values differ in length");
    std::vector<double> out(x.size());
    const auto last = xp.size() - 1;
    for (std::size_t k = 0; k < x.size(); ++k) {
        const double value = x[k];
        if (nanAfterEnd && value > xp[last]) { out[k] = std::numeric_limits<double>::quiet_NaN(); continue; }
        if (value < xp[0]) { out[k] = fp[0]; continue; }
        if (value > xp[last]) { out[k] = fp[last]; continue; }
        const auto j = static_cast<std::size_t>(std::upper_bound(xp.begin(), xp.end(), value) - xp.begin()) - 1;
        if (j == last || xp[j] == value) { out[k] = fp[j]; continue; }
        const double slope = (fp[j + 1] - fp[j]) / (xp[j + 1] - xp[j]);
        double result = slope * (value - xp[j]) + fp[j];
        if (std::isnan(result)) {
            result = slope * (value - xp[j + 1]) + fp[j + 1];
            if (std::isnan(result) && fp[j] == fp[j + 1]) result = fp[j];
        }
        out[k] = result;
    }
    return out;
}

std::vector<double> column(const c::NpyArray& array, std::size_t index) {
    const auto values = array.doubles();
    if (array.shape.size() == 1) {
        if (index) throw std::invalid_argument("figure_dataset: column of a vector");
        return values;
    }
    if (array.shape.size() != 2 || index >= array.shape[1]) throw std::invalid_argument("figure_dataset: column out of range");
    const auto rows = array.shape[0], columns = array.shape[1];
    std::vector<double> out(rows);
    for (std::size_t r = 0; r < rows; ++r) out[r] = values[r * columns + index];
    return out;
}

}  // namespace

std::filesystem::path runFigureDataset(const std::vector<std::filesystem::path>& sources, const json::Value& settings,
                                       const std::filesystem::path& root, const std::string& label) {
    std::map<std::string, c::NpzRecord> inputs;
    for (const auto& [name, path] : settings.at("inputs").object())
        inputs.emplace(name, c::readNpz(std::filesystem::u8path(path.string())));
    if (inputs.empty()) {
        for (std::size_t k = 0; k < sources.size(); ++k) inputs.emplace(std::to_string(k), c::readNpz(sources[k]));
    }
    const auto record = [&](const json::Value& spec) -> const c::NpzRecord& {
        const auto from = spec.at("from").string();
        if (!inputs.count(from)) throw std::invalid_argument("figure_dataset: unknown input '" + from + "'");
        return inputs.at(from);
    };
    c::NpzRecord out;
    for (const auto& spec : settings.at("arrays").array()) {
        const auto name = spec.at("name").string();
        if (name.empty()) throw std::invalid_argument("figure_dataset: every array needs a name");
        const auto& source = record(spec).at(spec.at("array").string());
        const bool matrix = !spec.at("columns").isNull();
        std::vector<std::size_t> indices;
        if (matrix) for (const auto& item : spec.at("columns").array()) indices.push_back(static_cast<std::size_t>(item.number()));
        else indices.push_back(static_cast<std::size_t>(spec.at("column").number(0.0)));
        const bool resample = !spec.at("onto").isNull();
        std::vector<double> clock, sourceClock;
        if (resample) {
            clock = record(spec.at("onto")).at(spec.at("onto").at("array").string()).doubles();
            sourceClock = record(spec).at(spec.at("clock").string("time")).doubles();
        }
        std::vector<std::vector<double>> series;
        for (const auto index : indices) {
            auto values = column(source, index);
            if (!spec.at("scale_before").isNull()) for (auto& value : values) value *= spec.at("scale_before").number();
            auto sampled = resample ? interpolate(clock, sourceClock, values, spec.at("nan_after_end").boolean(false))
                                    : std::move(values);
            if (!spec.at("scale").isNull()) for (auto& value : sampled) value *= spec.at("scale").number();
            if (!spec.at("offset").isNull()) for (auto& value : sampled) value += spec.at("offset").number();
            if (!spec.at("keep_times").isNull()) {  // NaN except at the samples nearest the listed times
                const auto times = resample ? clock : record(spec).at(spec.at("clock").string("time")).doubles();
                if (times.size() != sampled.size() || times.empty())
                    throw std::invalid_argument("figure_dataset: keep_times needs a clock matching the array");
                std::vector<double> kept(sampled.size(), std::numeric_limits<double>::quiet_NaN());
                for (const auto& item : spec.at("keep_times").array()) {
                    const double wanted = item.number();
                    auto at = static_cast<std::size_t>(std::lower_bound(times.begin(), times.end(), wanted) - times.begin());
                    if (at == times.size() || (at > 0 && wanted - times[at - 1] < times[at] - wanted)) --at;
                    kept[at] = sampled[at];
                }
                sampled = std::move(kept);
            }
            series.push_back(std::move(sampled));
        }
        const auto rows = series.front().size();
        if (!matrix) { out.arrays[name] = c::NpyArray::fromDoubles({rows}, series.front()); continue; }
        std::vector<double> joined(rows * series.size());
        for (std::size_t k = 0; k < series.size(); ++k)
            for (std::size_t r = 0; r < rows; ++r) joined[r * series.size() + k] = series[k][r];
        out.arrays[name] = c::NpyArray::fromDoubles({rows, series.size()}, joined);
    }
    if (out.arrays.empty()) throw std::invalid_argument("figure_dataset lists no arrays");
    out.arrays["metadata_json"] = c::NpyArray::fromString(json::dump(json::Value(json::Value::Object{
        {"format", json::Value("dynsyslab.derived/v1")}, {"method", json::Value("figure-dataset")},
        {"label", json::Value(label)}, {"dataset_settings", settings}})));
    const auto directory = objectFolder(root, "figure-dataset");
    std::filesystem::create_directories(directory);
    const auto path = directory / c::wide(safeName(label + "__" + timestamp()) + ".npz");
    c::writeNpzAtomicNew(path, out);
    return path;
}

}  // namespace dynsyslab::workflow::postprocess_detail
