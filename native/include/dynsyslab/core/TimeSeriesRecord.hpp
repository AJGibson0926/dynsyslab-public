#pragma once
// A time series as a pipeline data object: the samples (core::TimeSeries), what each channel is,
// and where the record came from. Methods that consume time series (e.g. H-DMD) see only this.
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/TimeSeries.hpp>

#include <optional>
#include <string>
#include <vector>

namespace dynsyslab::core {

// Pointer to a stored data object. The SHA-256 is the identity; the path is only a locator.
struct DataReference {
    std::string path;
    std::string sha256;
    std::string role;  // e.g. "time_series", "state_source"
};

// Every derived object records the method that produced it, its runtime settings and its inputs.
struct DataOrigin {
    std::string method;
    json::Value settings;
    std::vector<DataReference> inputs;
};

// One channel of a record. `observable` names the catalog observable that produced the channel
// ("state" components for a solver record, e.g. "echo" after a coordinate change); empty when the
// channel is a raw measurement of unknown provenance.
struct ChannelDescriptor {
    std::string name;
    std::string observable;
    std::string units;
};

struct TimeSeriesRecord {
    TimeSeries series;
    std::vector<ChannelDescriptor> channels;
    DataOrigin origin;
    json::Value system;  // configuration of the dynamical system that produced it, if any
    // Solver record whose states these channels were computed from (the record itself for a
    // solve). Absent for measured data: then no state domain exists.
    std::optional<DataReference> stateSource;
};

// Channel identities, in order; used to recognize interpolation domains.
[[nodiscard]] std::vector<std::string> channelObservables(const TimeSeriesRecord& record);

// Joins records sampled on one clock into a single multichannel record. Clocks must match exactly;
// nothing is resampled. Origin inputs and channels are concatenated in order.
[[nodiscard]] TimeSeriesRecord joinTimeSeriesRecords(const std::vector<TimeSeriesRecord>& records);

[[nodiscard]] json::Value dataOriginJson(const DataOrigin& origin);
[[nodiscard]] DataOrigin dataOriginFromJson(const json::Value& value);
[[nodiscard]] json::Value channelDescriptorsJson(const std::vector<ChannelDescriptor>& channels);
[[nodiscard]] std::vector<ChannelDescriptor> channelDescriptorsFromJson(const json::Value& value);

}  // namespace dynsyslab::core
