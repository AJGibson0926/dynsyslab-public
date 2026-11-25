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

c::TimeSeriesRecord timeSeriesRecord(const SourceData& source) {
    const auto meta = source.metadataJson.empty() ? J{} : json::parse(source.metadataJson);
    const auto method = meta.at("method").string();
    const auto kind = meta.at("coordinate_transform_kind").string();
    const bool hankelBasis = kind == "hankel_basis";
    const bool measured = meta.at("coordinate_provider").string() == "measured"
                       || (meta.at("value_space").string() == "history-coordinates" && !hankelBasis);
    const bool solve = method.empty() || method == "solve";
    const auto names = timeSeriesColumnNames(source);
    std::vector<c::ChannelDescriptor> channels;
    for (std::size_t k = 0; k < source.dimension; ++k) {
        std::string observable;
        if (solve || kind == "identity") observable = "state";
        else if (hankelBasis) observable = names[k];
        else if (kind == "observable") {
            observable = meta.at("coordinate_observable").string();
            if (source.dimension > 1 && !observable.empty()) observable += "[" + std::to_string(k) + "]";
        } else if (kind == "channels" && k < meta.at("coordinate_channels").array().size()) {
            const auto& channel = meta.at("coordinate_channels").array()[k];
            observable = channel.at("base").string();
            if (channel.at("derivative_source").string() != "none" && !channel.at("derivative_source").isNull())
                observable += "'";
        }
        if (measured) observable.clear();
        channels.push_back({names[k], observable, ""});
    }
    c::TimeSeriesRecord record{c::TimeSeries(source.time, source.values, source.dimension), std::move(channels),
                               {}, meta.at("system_configuration"), std::nullopt};
    record.origin.method = solve ? "solve" : method;
    record.origin.inputs.push_back(dataReference(source.path, "time_series"));
    // A state domain exists whenever a producing solve does, even if a channel was computed
    // numerically ("measured" provider) rather than as an analytic observable.
    if (solve && !measured) record.stateSource = record.origin.inputs.front();
    else if (!meta.at("state_source_archive").string().empty())
        record.stateSource = dataReference(std::filesystem::u8path(meta.at("state_source_archive").string()), "state_source");
    else if (!solve && !hankelBasis && !meta.at("source_archive").string().empty())
        record.stateSource = dataReference(std::filesystem::u8path(meta.at("source_archive").string()), "state_source");
    return record;
}

c::TimeSeriesRecord loadTimeSeriesRecords(const std::vector<std::filesystem::path>& paths) {
    std::vector<c::TimeSeriesRecord> records;
    for (const auto& path : paths) records.push_back(timeSeriesRecord(loadSource(path)));
    return c::joinTimeSeriesRecords(records);
}

SourceData sourceView(const c::TimeSeriesRecord& record, const SourceData& like) {
    const auto& series = record.series;
    SourceData view;
    view.path = like.path;
    view.samples = series.sampleCount();
    view.dimension = series.componentCount();
    view.time.assign(series.timeData(), series.timeData() + view.samples);
    view.values.assign(series.valuesData(), series.valuesData() + view.samples * view.dimension);
    auto meta = like.metadataJson.empty() ? J::Object{} : json::parse(like.metadataJson).object();
    J::Array names;
    for (const auto& channel : record.channels) names.emplace_back(channel.name);
    meta["coordinate_names"] = names;
    view.metadataJson = json::dump(J(meta));
    return view;
}

}  // namespace dynsyslab::workflow::postprocess_detail
