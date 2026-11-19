#pragma once

#include <dynsyslab/core/Archive.hpp>

#include <optional>

namespace dsl {

// A time window in the archive's stored clock units (the values of its `time` array).
struct TimeWindow {
    std::optional<double> start, stop;
    bool active() const { return start.has_value() || stop.has_value(); }
};

// Restricts every array aligned with a recognized clock (time, controlled_time,
// control_time, free_time, estimator_time) to the rows whose clock value lies in
// [start, stop]. An array is aligned with a clock when its leading dimension equals
// that clock's length; other arrays are kept whole. Every Viewer route (time series,
// phase space, Fourier spectra, estimator panels) built from the result therefore
// respects the window.
dynsyslab::core::NpzRecord windowRecord(const dynsyslab::core::NpzRecord& record,
                                        const TimeWindow& window);

}  // namespace dsl
