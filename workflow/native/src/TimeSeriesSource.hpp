#pragma once
// Adapter from archived records to core::TimeSeriesRecord, the input type of time-series methods.
#include "PostprocessInternal.hpp"
#include <dynsyslab/core/TimeSeriesRecord.hpp>

namespace dynsyslab::workflow::postprocess_detail {

// Wraps a loaded (unit-normalized) source: samples, channel identities and origin. A solver record
// is its own state source; a coordinate change points to the solver record it was computed from.
[[nodiscard]] dynsyslab::core::TimeSeriesRecord timeSeriesRecord(const SourceData& source);

// Loads every path and joins them on their common clock.
[[nodiscard]] dynsyslab::core::TimeSeriesRecord loadTimeSeriesRecords(const std::vector<std::filesystem::path>& paths);

// A SourceData view of a record, for the shared Hankel window builder. Carries channel names only:
// no system state or other side data.
[[nodiscard]] SourceData sourceView(const dynsyslab::core::TimeSeriesRecord& record, const SourceData& like);

[[nodiscard]] dynsyslab::core::DataReference dataReference(const std::filesystem::path& path, const std::string& role);

}  // namespace dynsyslab::workflow::postprocess_detail
