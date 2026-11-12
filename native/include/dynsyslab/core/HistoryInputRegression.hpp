#pragma once

#include <array>
#include <complex>
#include <cstddef>
#include <vector>
#include <dynsyslab/core/KoopmanModel.hpp>

namespace dynsyslab::core {

// Measurement-only boundary. Coordinates and psi are obtained from frozen
// echo histories; no system, state, derivatives or model callbacks are accepted.
struct HistoryInputSample {
    std::array<double, 2> coordinates{};
    std::complex<double> psi{}, nextPsi{};
    double input{};
};

struct HistoryInputMap {
    std::size_t degree{}; // 0 constant, 1 conjugate-complete bilinear in psi
    std::array<double, 2> coordinateScale{1.0, 1.0};
    std::vector<std::complex<double>> coefficients;
    double sampleTime{};
    std::complex<double> multiplier{};
    std::size_t samples{}, rank{};
    double conditionNumber{};
    [[nodiscard]] std::complex<double> evaluate(std::complex<double> psi) const;
    [[nodiscard]] ApproximateControlMap canonicalRepresentation() const;
};

struct HistoryInputMetrics {
    double relativeResidual{}, rmse{}, unforcedBaselineRmse{};
};

// z[k+1] - exp(lambda*dt) z[k] = B_d(xi[k]) u[k].
// This is a sampled response fit, NOT an exact instantaneous Jacobian.
[[nodiscard]] HistoryInputMap fitHistoryInputMap(
    const std::vector<HistoryInputSample>& samples, std::complex<double> multiplier,
    double sampleTime, std::size_t degree, std::array<double, 2> coordinateScale,
    double relativeRidge = 1e-8);
[[nodiscard]] HistoryInputMetrics scoreHistoryInputMap(
    const HistoryInputMap& map, const std::vector<HistoryInputSample>& samples);

} // namespace dynsyslab::core
