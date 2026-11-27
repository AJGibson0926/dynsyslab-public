#include <dynsyslab/core/PathSampling.hpp>

#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>

namespace dynsyslab::core {

std::vector<double> resampleByArcLength(const std::vector<double>& states, std::size_t samples, std::size_t n,
                                        std::size_t count) {
    if (states.size() != samples * n || samples < 2 || count < 2)
        throw std::invalid_argument("arc-length resampling needs a path of two or more states");
    std::vector<double> arc(samples, 0.0);
    for (std::size_t k = 1; k < samples; ++k) {
        double d = 0.0;
        for (std::size_t i = 0; i < n; ++i) d += std::pow(states[k * n + i] - states[(k - 1) * n + i], 2);
        arc[k] = arc[k - 1] + std::sqrt(d);
    }
    if (!(arc.back() > 0.0)) throw std::invalid_argument("arc-length resampling: the path does not move");
    std::vector<double> result;
    result.reserve(count * n);
    std::size_t segment = 1;
    for (std::size_t m = 0; m < count; ++m) {
        const double s = arc.back() * double(m) / double(count - 1);
        while (segment + 1 < samples && arc[segment] < s) ++segment;
        const double width = arc[segment] - arc[segment - 1];
        const double w = width > 0.0 ? std::clamp((s - arc[segment - 1]) / width, 0.0, 1.0) : 0.0;
        for (std::size_t i = 0; i < n; ++i)
            result.push_back((1 - w) * states[(segment - 1) * n + i] + w * states[segment * n + i]);
    }
    return result;
}

std::vector<double> tubeStates(const std::vector<double>& path, std::size_t samples, std::size_t n, double radius,
                               std::size_t perPoint, std::uint64_t seed) {
    if (path.size() != samples * n || !samples) throw std::invalid_argument("tube sampling needs a path");
    if (!(radius >= 0.0)) throw std::invalid_argument("tube radius must be nonnegative");
    std::vector<double> scale(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        double low = path[i], high = path[i];
        for (std::size_t k = 0; k < samples; ++k) { low = std::min(low, path[k * n + i]); high = std::max(high, path[k * n + i]); }
        scale[i] = std::max(high - low, 1e-12);
    }
    std::mt19937_64 random(seed);
    std::normal_distribution<double> normal;
    std::vector<double> result = path;
    for (std::size_t k = 0; k < samples; ++k)
        for (std::size_t p = 0; p < perPoint; ++p) {
            std::vector<double> direction(n);
            double norm = 0.0;
            for (auto& d : direction) { d = normal(random); norm += d * d; }
            norm = std::sqrt(norm);
            for (std::size_t i = 0; i < n; ++i)
                result.push_back(path[k * n + i] + radius * scale[i] * (norm > 0 ? direction[i] / norm : 0.0));
        }
    return result;
}

}  // namespace dynsyslab::core
