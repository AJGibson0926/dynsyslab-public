#include <dynsyslab/core/TimeSeriesRecord.hpp>

#include <stdexcept>

namespace dynsyslab::core {

std::vector<std::string> channelObservables(const TimeSeriesRecord& record) {
    std::vector<std::string> names;
    for (const auto& channel : record.channels) names.push_back(channel.observable);
    return names;
}

TimeSeriesRecord joinTimeSeriesRecords(const std::vector<TimeSeriesRecord>& records) {
    if (records.empty()) throw std::invalid_argument("joinTimeSeriesRecords needs at least one record");
    if (records.size() == 1) return records.front();
    const auto& first = records.front().series;
    const std::size_t samples = first.sampleCount();
    std::size_t width = 0;
    for (const auto& record : records) {
        const auto& series = record.series;
        if (series.sampleCount() != samples)
            throw std::invalid_argument("joined time series must share one clock (sample counts differ)");
        for (std::size_t k = 0; k < samples; ++k)
            if (series.time(k) != first.time(k))
                throw std::invalid_argument("joined time series must share one clock (times differ); "
                                            "resample them first");
        if (record.inputs)
            throw std::invalid_argument("joined time series must not carry inputs; drop them before joining");
        if (record.channels.size() != series.componentCount())
            throw std::invalid_argument("time-series record channels do not match its components");
        width += series.componentCount();
    }
    std::vector<double> time(first.timeData(), first.timeData() + samples), values(samples * width);
    TimeSeriesRecord joined;
    joined.origin.method = "join";
    std::size_t offset = 0;
    for (const auto& record : records) {
        const auto& series = record.series;
        for (std::size_t k = 0; k < samples; ++k)
            for (std::size_t c = 0; c < series.componentCount(); ++c)
                values[k * width + offset + c] = series.value(k, c);
        offset += series.componentCount();
        joined.channels.insert(joined.channels.end(), record.channels.begin(), record.channels.end());
        joined.origin.inputs.insert(joined.origin.inputs.end(), record.origin.inputs.begin(),
                                    record.origin.inputs.end());
        if (joined.system.isNull()) joined.system = record.system;
        if (!joined.stateSource && record.stateSource) joined.stateSource = record.stateSource;
        else if (joined.stateSource && record.stateSource && joined.stateSource->sha256 != record.stateSource->sha256)
            joined.stateSource.reset();  // channels come from different trajectories: no single state domain
    }
    joined.series = TimeSeries(std::move(time), std::move(values), width);
    return joined;
}

json::Value dataOriginJson(const DataOrigin& origin) {
    json::Value::Array inputs;
    for (const auto& input : origin.inputs)
        inputs.emplace_back(json::Value::Object{{"path", input.path}, {"sha256", input.sha256}, {"role", input.role}});
    return json::Value::Object{{"schema", "dynsyslab.data-origin/v1"}, {"method", origin.method},
                               {"settings", origin.settings}, {"inputs", inputs}};
}

DataOrigin dataOriginFromJson(const json::Value& value) {
    DataOrigin origin{value.at("method").string(), value.at("settings"), {}};
    if (!value.at("inputs").isNull())
        for (const auto& input : value.at("inputs").array())
            origin.inputs.push_back({input.at("path").string(), input.at("sha256").string(), input.at("role").string()});
    return origin;
}

json::Value channelDescriptorsJson(const std::vector<ChannelDescriptor>& channels) {
    json::Value::Array out;
    for (const auto& channel : channels)
        out.emplace_back(json::Value::Object{{"name", channel.name}, {"observable", channel.observable},
                                             {"units", channel.units}});
    return out;
}

std::vector<ChannelDescriptor> channelDescriptorsFromJson(const json::Value& value) {
    std::vector<ChannelDescriptor> channels;
    if (value.isNull()) return channels;
    for (const auto& item : value.array())
        channels.push_back({item.at("name").string(), item.at("observable").string(), item.at("units").string()});
    return channels;
}

}  // namespace dynsyslab::core
