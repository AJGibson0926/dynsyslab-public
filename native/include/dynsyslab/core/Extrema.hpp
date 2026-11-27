#pragma once
// Sample extrema used to summarize long records.
#include <cstddef>
#include <utility>
#include <vector>

namespace dynsyslab::core {

// Indices of the local maxima of `values` in its last `tailFraction` of samples (0 < f <= 1):
// i with v[i-1] <= v[i] >= v[i+1], after the transient start = max(1, floor((1 - f) n)). These
// are the points of an orbit (peak-amplitude) bifurcation diagram: one point per oscillation peak,
// so periodic orbits give a few points and chaotic ones a band (Strogatz, Nonlinear Dynamics and
// Chaos, 2nd ed., 2015, sec. 10.2; Parlitz & Lauterborn, Phys. Rev. A 36, 1428, 1987).
[[nodiscard]] std::vector<std::size_t> postTransientMaxima(const std::vector<double>& values,
                                                           double tailFraction = 0.25);

// A frequency range [0, upper] that holds a spectrum's dominant content: upper is the larger of
// the frequency below which `mass` (default 99.5%) of the non-negative amplitude lies and 1.25
// times the peak frequency. Only positive frequencies count. Empty when nothing is positive.
[[nodiscard]] std::vector<double> dominantFrequencyRange(const std::vector<double>& frequencies,
                                                         const std::vector<double>& amplitudes,
                                                         double mass = 0.995);

}  // namespace dynsyslab::core
