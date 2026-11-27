#include "SolveArchive.hpp"

#include <dynsyslab/core/TimeSeriesOperations.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace dynsyslab::workflow::detail {

void truncateOutputRows(core::NpzRecord& archive, std::size_t sampleCount, std::size_t rows) {
    for (auto& [name, array] : archive.arrays) {
        if (array.shape.empty() || array.shape[0] != sampleCount) continue;
        const std::size_t rowBytes = array.data.size() / array.shape[0];
        array.data.resize(rows * rowBytes);
        array.shape[0] = rows;
    }
}

json::Value resampleOutputArrays(core::NpzRecord& archive, const std::vector<double>& time,
                                 std::size_t rows, const ResampleSettings& settings) {
    const std::size_t samples = time.size();
    if (rows < 1) throw std::invalid_argument("resampling needs at least one computed sample");
    const std::vector<double> source(time.begin(), time.begin() + static_cast<std::ptrdiff_t>(rows));
    std::vector<std::string> names;
    for (const auto& [name, array] : archive.arrays)
        if (name != "time" && array.type == core::ArrayType::float64 && !array.shape.empty()
            && array.shape[0] == samples && name.rfind("raw_", 0) != 0)
            names.push_back(name);
    core::ResampleOptions options{settings.method, core::Extrapolation::nan};
    for (const auto& name : names) {
        const auto original = archive.arrays.at(name);
        const auto values = original.doubles();
        const std::size_t width = samples ? values.size() / samples : 0;
        std::vector<double> kept(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(rows * width));
        const auto resampled = core::TimeSeries(source, std::move(kept), width).resample(settings.times, options);
        auto shape = original.shape;
        shape[0] = settings.times.size();
        archive.arrays["raw_" + name] = original;
        archive.arrays[name] = core::NpyArray::fromDoubles(
            std::move(shape), resampled.valuesData(), resampled.sampleCount() * width);
    }
    archive.arrays["raw_time"] = core::NpyArray::fromDoubles({samples}, time);
    archive.arrays["time"] = core::NpyArray::fromDoubles({settings.times.size()}, settings.times);
    return json::Value::Object{{"method", std::string(core::interpolationName(settings.method))},
                               {"raw_samples", static_cast<double>(samples)},
                               {"output_samples", static_cast<double>(settings.times.size())},
                               {"raw_arrays_prefix", std::string("raw_")}};
}

std::filesystem::path partialArchivePath(const std::filesystem::path& archive) {
    auto name = archive.filename().wstring();
    if (name.size() > 4 && name.substr(name.size() - 4) == L".npz") name.resize(name.size() - 4);
    return archive.parent_path() / (name + L".partial.npz");
}

void writePartialArchive(const std::filesystem::path& path, const core::IntegrationCheckpoint& point,
                         const std::string& metadataJson) {
    const std::size_t rows = point.sample + 1;
    core::NpzRecord record;
    record.arrays["time"] = core::NpyArray::fromDoubles({rows}, point.times, rows);
    record.arrays["values"] = core::NpyArray::fromDoubles(
        {rows, point.stateDimension}, point.values, rows * point.stateDimension);
    if (point.controls)
        record.arrays["control_values"] = core::NpyArray::fromDoubles(
            {rows, point.inputDimension}, point.controls, rows * point.inputDimension);
    record.arrays["metadata_json"] = core::NpyArray::fromString(metadataJson);
    core::writeNpzAtomic(path, record);
}

void writePartialSeries(const std::filesystem::path& path, const double* time, std::size_t rows,
                        const double* values, std::size_t width, const std::string& metadataJson,
                        const std::map<std::string, core::NpyArray>& extra) {
    core::NpzRecord record;
    record.arrays["time"] = core::NpyArray::fromDoubles({rows}, time, rows);
    record.arrays["values"] = core::NpyArray::fromDoubles({rows, width}, values, rows * width);
    for (const auto& [name, array] : extra) record.arrays[name] = array;
    record.arrays["metadata_json"] = core::NpyArray::fromString(metadataJson);
    core::writeNpzAtomic(path, record);
}

double checkpointInterval(const json::Value& job) {
    const double interval = job.at("history").at("checkpoint_interval").number(0.0);
    if (!(interval >= 0.0) || !std::isfinite(interval))
        throw std::invalid_argument("history.checkpoint_interval must be nonnegative");
    return interval;
}

}  // namespace dynsyslab::workflow::detail
