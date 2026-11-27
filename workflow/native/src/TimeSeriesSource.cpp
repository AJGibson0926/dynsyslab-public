#include "TimeSeriesSource.hpp"
#include "CoordinateChannels.hpp"
#include <dynsyslab/core/FileDigest.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace c = dynsyslab::core;
using J = json::Value;

c::DataReference dataReference(const std::filesystem::path& path, const std::string& role) {
    const auto absolute = std::filesystem::absolute(path).lexically_normal();
    return {absolute.generic_u8string(), c::sha256File(absolute), role};
}

// A referenced solver record that cannot be opened keeps its locator with an empty digest; a method
// that needs its states fails when it opens it.
c::DataReference stateSourceReference(const std::string& path) {
    const auto file = std::filesystem::u8path(path);
    std::error_code error;
    if (!std::filesystem::is_regular_file(file, error))
        return {std::filesystem::absolute(file).lexically_normal().generic_u8string(), "", "state_source"};
    return dataReference(file, "state_source");
}

void describeRecord(SourceData& source) {
    const auto meta = source.metadataJson.empty() ? J{} : json::parse(source.metadataJson);
    const auto method = meta.at("method").string();
    const auto kind = meta.at("coordinate_transform_kind").string();
    const bool hankelBasis = kind == "hankel_basis";
    const bool measured = meta.at("coordinate_provider").string() == "measured"
                       || (meta.at("value_space").string() == "history-coordinates" && !hankelBasis);
    const bool solve = method.empty() || method == "solve";
    const auto names = timeSeriesColumnNames(source);
    std::vector<c::ChannelDescriptor> channels;
    for (std::size_t k = 0; k < source.dimension(); ++k) {
        std::string observable;
        if (solve || kind == "identity") observable = "state";
        else if (hankelBasis) observable = names[k];
        else if (kind == "observable") {
            observable = meta.at("coordinate_observable").string();
            if (source.dimension() > 1 && !observable.empty()) observable += "[" + std::to_string(k) + "]";
        } else if (kind == "channels" && k < meta.at("coordinate_channels").array().size()) {
            const auto& channel = meta.at("coordinate_channels").array()[k];
            observable = channel.at("base").string();
            if (channel.at("derivative_source").string() != "none" && !channel.at("derivative_source").isNull())
                observable += "'";
        }
        if (measured) observable.clear();
        channels.push_back({names[k], observable, ""});
    }
    auto& record = source.record;
    record.channels = std::move(channels);
    record.origin = {};
    record.system = meta.at("system_configuration");
    record.stateSource.reset();
    record.origin.method = solve ? "solve" : method;
    record.origin.inputs.push_back(dataReference(source.path, "time_series"));
    // A state domain exists whenever a producing solve does, even if a channel was computed
    // numerically ("measured" provider) rather than as an analytic observable.
    if (solve && !measured) record.stateSource = record.origin.inputs.front();
    else if (!meta.at("state_source_archive").string().empty())
        record.stateSource = stateSourceReference(meta.at("state_source_archive").string());
    else if (!solve && !hankelBasis && !meta.at("source_archive").string().empty())
        record.stateSource = stateSourceReference(meta.at("source_archive").string());
}

c::TimeSeriesRecord loadTimeSeriesRecords(const std::vector<std::filesystem::path>& paths) {
    std::vector<c::TimeSeriesRecord> records;
    for (const auto& path : paths) {
        auto record = loadSource(path).record;
        record.inputs.reset(); record.inputsOnIntervals = false;  // the samples only
        records.push_back(std::move(record));
    }
    return c::joinTimeSeriesRecords(records);
}

SourceData joinedTimeSeries(const SourceData& full, const J& settings) {
    // Time series only: inputs and other side data are not part of the joined record.
    auto part = [](c::TimeSeriesRecord record) { record.inputs.reset(); record.inputsOnIntervals = false; return record; };
    std::vector<c::TimeSeriesRecord> parts{part(full.record)};
    for (const auto& path : settings.at("additional_time_series").array())
        parts.push_back(part(loadSource(std::filesystem::u8path(path.string())).record));
    SourceData view{full.path, {}, c::joinTimeSeriesRecords(parts)};
    auto meta = full.metadataJson.empty() ? J::Object{} : json::parse(full.metadataJson).object();
    J::Array names;
    for (const auto& channel : view.record.channels) names.emplace_back(channel.name);
    meta["coordinate_names"] = names;
    view.metadataJson = json::dump(J(meta));
    return view;
}

}  // namespace dynsyslab::workflow::postprocess_detail
