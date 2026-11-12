#pragma once
#include <cmath>
namespace dynsyslab::core::detail {
struct HodgkinHuxleyRates {
    double am{}, bm{}, ah{}, bh{}, an{}, bn{};
};
inline HodgkinHuxleyRates hodgkinHuxleyRates(double voltage) {
    const auto ratio = [](double x, double scale, double coefficient) {
        return std::abs(x) < 1e-8 ? coefficient * scale : coefficient * x / -std::expm1(-x / scale);
    };
    return {ratio(voltage + 40, 10, .1),          4 * std::exp(-(voltage + 65) / 18),
            .07 * std::exp(-(voltage + 65) / 20), 1 / (1 + std::exp(-(voltage + 35) / 10)),
            ratio(voltage + 55, 10, .01),         .125 * std::exp(-(voltage + 65) / 80)};
}
} // namespace dynsyslab::core::detail
