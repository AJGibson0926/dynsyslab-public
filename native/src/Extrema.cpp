// Record extrema (core/Extrema.hpp).
#include <dynsyslab/core/Extrema.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::core {

std::vector<std::size_t> postTransientMaxima(const std::vector<double>& values, double tailFraction) {
    if (!(tailFraction > 0.0 && tailFraction <= 1.0))
        throw std::invalid_argument("tail fraction must lie in (0, 1]");
    const std::size_t n = values.size();
    std::vector<std::size_t> maxima;
    if (n < 3) return maxima;
    const auto start = std::max<std::size_t>(1, static_cast<std::size_t>(std::floor((1.0 - tailFraction) * n)));
    for (std::size_t i = start + 1; i + 1 < n; ++i)
        if (values[i] >= values[i - 1] && values[i] >= values[i + 1] && std::isfinite(values[i]))
            maxima.push_back(i);
    return maxima;
}

std::vector<double> dominantFrequencyRange(const std::vector<double>& frequencies,
                                           const std::vector<double>& amplitudes, double mass) {
    if (frequencies.size() != amplitudes.size())
        throw std::invalid_argument("frequencies and amplitudes differ in length");
    if (!(mass > 0.0 && mass <= 1.0)) throw std::invalid_argument("spectral mass must lie in (0, 1]");
    double total = 0.0, peakValue = -1.0, peak = 0.0, last = 0.0;
    for (std::size_t i = 0; i < frequencies.size(); ++i) {
        if (!(frequencies[i] > 0.0) || !std::isfinite(amplitudes[i])) continue;
        const double weight = std::max(amplitudes[i], 0.0);
        total += weight; last = frequencies[i];
        if (weight > peakValue) { peakValue = weight; peak = frequencies[i]; }
    }
    if (last <= 0.0) return {};
    if (total <= 0.0) return {0.0, last};
    double running = 0.0, cutoff = last;
    for (std::size_t i = 0; i < frequencies.size(); ++i) {
        if (!(frequencies[i] > 0.0) || !std::isfinite(amplitudes[i])) continue;
        running += std::max(amplitudes[i], 0.0);
        if (running >= mass * total) { cutoff = frequencies[i]; break; }
    }
    return {0.0, std::min(last, std::max(cutoff, 1.25 * peak))};
}

}  // namespace dynsyslab::core
