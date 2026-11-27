// Post-transient maxima, dominant spectral band and observability mode participation.
#include <dynsyslab/core/Extrema.hpp>
#include <dynsyslab/core/Observability.hpp>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
void require(bool condition, const char* message) { if (!condition) throw std::runtime_error(message); }
constexpr double kPi = 3.14159265358979323846;
}

int main() {
    try {
        using namespace dynsyslab::core;
        // A decaying transient then a unit sine of period 100 samples: one maximum per period in the tail.
        std::vector<double> values(1000);
        for (std::size_t i = 0; i < values.size(); ++i)
            values[i] = std::sin(2 * kPi * i / 100.0) + (i < 200 ? 5.0 * std::exp(-0.02 * i) : 0.0);
        const auto maxima = postTransientMaxima(values, 0.5);
        require(maxima.size() == 5, "five maxima in the last half");
        for (const auto i : maxima) require(i >= 500 && std::abs(values[i] - 1.0) < 1e-3, "maxima at the crests");
        require(postTransientMaxima({1, 2}, 1.0).empty(), "short records have no maxima");
        bool rejected = false;
        try { (void)postTransientMaxima(values, 0.0); } catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "a zero tail fraction is rejected");

        // A spectrum peaked at 2 with all mass below 3: the band ends at max(cutoff, 1.25 peak).
        std::vector<double> f, a;
        for (int k = 0; k <= 100; ++k) { f.push_back(0.1 * k); a.push_back(f.back() > 0 && f.back() <= 3 ? std::exp(-std::pow(f.back() - 2, 2)) : 0.0); }
        const auto band = dominantFrequencyRange(f, a);
        require(band.size() == 2 && band[0] == 0.0 && band[1] >= 2.5 && band[1] <= 3.0 + 1e-12, "dominant band");
        require(dominantFrequencyRange({0.0}, {1.0}).empty(), "no positive frequency, no band");

        // An orthonormal rotation: participation of a full partition sums to one in every direction.
        ObservabilityResult result;
        result.stateDimension = 3;
        const double c = std::cos(0.3), s = std::sin(0.3);
        result.rightDirections = {c, -s, 0, s, c, 0, 0, 0, 1};
        const auto participation = observabilityModeParticipation(result, {{0, 1}, {2}});
        require(participation.size() == 6, "groups x directions");
        for (std::size_t j = 0; j < 3; ++j)
            require(std::abs(participation[j] + participation[3 + j] - 1.0) < 1e-12, "columns sum to one");
        require(std::abs(participation[2]) < 1e-12 && std::abs(participation[5] - 1.0) < 1e-12, "third direction is mode 2");
        rejected = false;
        try { (void)observabilityModeParticipation(result, {{0, 1}, {1}}); } catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "overlapping groups are rejected");
        std::cout << "dynsyslab_extrema_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_extrema_tests: FAIL: " << error.what() << '\n';
        return 1;
    }
}
