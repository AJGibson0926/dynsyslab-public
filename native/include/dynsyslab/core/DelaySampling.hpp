#pragma once
// Delayed observations at physical-time (fractional) lags, for delay coordinates whose taps do not
// fall on the sample grid (Takens, Lect. Notes Math. 898:366, 1981; delay embedding with real-valued
// lags: Sauer, Yorke & Casdagli, J. Stat. Phys. 65:579, 1991).
//
//   delaySamples      y(t_e - tau_l) of every channel, linear between samples. Causal: a tap never
//                     reads a sample after its endpoint, and a window that would start before the
//                     record is refused. Lags within roundoff of a tick read that sample exactly.
//   sampleDelayJet    (y, y', y'')(t_e - tau_l) from a uniformly sampled jet (value, rate,
//                     acceleration) with the C^2 quintic Hermite interpolant on each interval (the
//                     unique quintic matching value and first two derivatives at both ends).
#include <dynsyslab/core/TimeSeries.hpp>

#include <vector>

namespace dynsyslab::core {

// endpoints x lags x components, row-major.
[[nodiscard]] std::vector<double> delaySamples(const TimeSeries& series, const std::vector<double>& endpoints,
                                               const std::vector<double>& lags);

// `jet` holds samples x 3 (value, rate, acceleration) on a uniform clock; endpoints x lags x 3.
[[nodiscard]] std::vector<double> sampleDelayJet(const TimeSeries& jet, const std::vector<double>& endpoints,
                                                 const std::vector<double>& lags);

}  // namespace dynsyslab::core
