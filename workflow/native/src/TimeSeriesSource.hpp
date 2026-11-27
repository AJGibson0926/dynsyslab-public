#pragma once
// Adapter from archived records to core::TimeSeriesRecord, the input type of time-series methods.
#include "PostprocessInternal.hpp"
#include <dynsyslab/core/TimeSeriesRecord.hpp>

namespace dynsyslab::workflow::postprocess_detail {

// Fills source.record's channel identities, origin, system and state source from its path and
// (unit-normalized) metadata; loadSource calls it. A solver record is its own state source; a
// coordinate change points to the solver record it was computed from.
void describeRecord(SourceData& source);

// Loads every path and joins them on their common clock.
[[nodiscard]] dynsyslab::core::TimeSeriesRecord loadTimeSeriesRecords(const std::vector<std::filesystem::path>& paths);

// `full` joined with settings.additional_time_series on their common clock, for the shared Hankel
// window builder: the samples and channel names only (no inputs, system state or other side data).
[[nodiscard]] SourceData joinedTimeSeries(const SourceData& full, const json::Value& settings);

[[nodiscard]] dynsyslab::core::DataReference dataReference(const std::filesystem::path& path, const std::string& role);

}  // namespace dynsyslab::workflow::postprocess_detail
