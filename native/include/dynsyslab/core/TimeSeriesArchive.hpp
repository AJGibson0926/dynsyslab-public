#pragma once

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/TimeSeries.hpp>

#include <string>
#include <vector>

namespace dynsyslab::core {

struct TimeSeriesArchiveSelection {
    std::vector<std::string> timeArrays{"time"};
    std::vector<std::string> valueArrays{"values"};
};

// Materializes one selected archive record as an owning TimeSeries. It never
// reads or writes derived Fourier data.
[[nodiscard]] TimeSeries loadTimeSeries(
    const NpzRecord& archive, const TimeSeriesArchiveSelection& selection = {});

// Writes only the authoritative time and sample arrays. An empty time-array
// name reuses a time axis already written to the record.
void storeTimeSeries(NpzRecord& archive, const TimeSeries& series,
                     const std::string& timeArray = "time",
                     const std::string& valueArray = "values");

}  // namespace dynsyslab::core
