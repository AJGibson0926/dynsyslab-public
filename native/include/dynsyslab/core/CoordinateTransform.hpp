#pragma once
#include <dynsyslab/core/TimeSeries.hpp>
#include <dynsyslab/core/Observable.hpp>
#include <dynsyslab/core/HistoryObservation.hpp>
#include <optional>

namespace dynsyslab::core {
struct CoordinateSeries {
    TimeSeries coordinates;
    std::optional<TimeSeries> derivatives;
};
// Stateless analytic observation. Rates, if supplied, are pushed through D h;
// the caller decides whether those rates were measured or evaluated from a model.
CoordinateSeries transformCoordinates(const TimeSeries& source, const Observable& map,
    const TimeSeries* sourceRates = nullptr);
// Stateful causal observation. No SVD/refit; outputs are timestamped at window END.
// Applying the same linear map to measurement rates differentiates the history
// coordinates along a trajectory. It does not prove autonomous closure in them.
CoordinateSeries transformHistoryCoordinates(const TimeSeries& scalarMeasurement,
    const FrozenHistoryTransform& transform, const TimeSeries* measurementRates = nullptr);
// Returns (measurement, signed trailing integral) at each complete window end.
// Duration is in the input clock's units. Interpolates linearly between samples.
TimeSeries windowIntegralCoordinates(const TimeSeries& scalarMeasurement, double duration);
// Fourth-order centered, uniform-clock derivative. Drops two endpoints on each
// side; never invents endpoint rates or silently differentiates across records.
TimeSeries differentiateTimeSeries(const TimeSeries& source);
TimeSeries trimTimeSeries(const TimeSeries& source, std::size_t first, std::size_t count);
void requireAlignedTimeSeries(const TimeSeries& left, const TimeSeries& right);
}
