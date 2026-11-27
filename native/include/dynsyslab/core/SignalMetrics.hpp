#pragma once
// Scalar-signal metrics for scoring controlled runs (WF control_score).
//   fitTone:        least-squares fit of offset + sum_k a_k cos(k w t) + b_k sin(k w t) at a known
//                   fundamental (the linear sine-fitting of IEEE Std 1057-2017, Sec. 4.1.3, with
//                   harmonics); amplitudes, phases (of sin(k w t + phi)), residual RMS and the
//                   harmonic content sqrt(sum_{k>=2} A_k^2) / A_1 (total harmonic distortion).
//   settlingTime:   the first time after which |e| stays within the band (NaN if never).
//   lowPass:        first-order low-pass filter, exact for a held input over each sample:
//                   y_k = y_{k-1} + (1 - exp(-2 pi f dt)) (x_k - y_{k-1}).
//   autocorrelation: sample autocorrelation at a lag (whiteness of innovations: Box & Jenkins,
//                   Time Series Analysis, 1970, Ch. 8; near zero for a white sequence).
#include <cstddef>
#include <vector>

namespace dynsyslab::core {

struct ToneFit {
    double offset{};
    std::vector<double> amplitudes;  // harmonic 1..H
    std::vector<double> phases;      // radians, as sin(k w t + phase)
    double residualRms{};
    double harmonicContent{};
};

[[nodiscard]] ToneFit fitTone(const std::vector<double>& time, const std::vector<double>& values,
                              double frequencyHz, std::size_t harmonics = 1);
[[nodiscard]] double settlingTime(const std::vector<double>& time, const std::vector<double>& error,
                                  double band);
[[nodiscard]] std::vector<double> lowPass(const std::vector<double>& time,
                                          const std::vector<double>& values, double cutoffHz);
[[nodiscard]] double autocorrelation(const std::vector<double>& values, std::size_t lag);

}  // namespace dynsyslab::core
