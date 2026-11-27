#pragma once
// Scores a forecast against a reference on its own clock (vector-field SINDy validation).
//
//   forecastScore   the prediction interpolated linearly onto the reference times inside its own
//                   span (never extrapolated): coverage of the requested interval, RMS and maximum
//                   errors (optionally trapezoid time-weighted), and, with a prominence, peak matching:
//                   reference peaks of one channel (local maxima whose topographic prominence reaches
//                   the threshold, as scipy.signal.find_peaks) each own the interval between the
//                   midpoints to its neighbours, and the highest forecast peak there is its match
//                   (time error, relative height). Peaks are found on each series' native clock.
//   kinematicMismatch  for coordinates holding a value and its rate, the model's chain rule d(value)/dt
//                   = f_value(z) against the rate coordinate: ||f_value - z_rate|| / ||z_rate||.
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/SindyDynamics.hpp>
#include <dynsyslab/core/TimeSeries.hpp>

#include <cstddef>
#include <optional>
#include <vector>

namespace dynsyslab::core {

struct ForecastScoreOptions {
    bool timeWeighted{false};
    std::optional<double> prominence;  // peak matching when set (>= 0)
    std::size_t channel{0};            // the channel peaks are matched on
};
[[nodiscard]] json::Value forecastScore(const TimeSeries& reference, const TimeSeries& prediction,
                                        const ForecastScoreOptions& options = {});

// Local maxima of `values` with topographic prominence at least `prominence` (sample indices).
[[nodiscard]] std::vector<std::size_t> prominentPeaks(const std::vector<double>& values, double prominence);

[[nodiscard]] double kinematicMismatch(const SindyDynamicsModel& model, const TimeSeries& coordinates,
                                       std::size_t valueChannel, std::size_t rateChannel);

}  // namespace dynsyslab::core
