#pragma once

#include <complex>
#include <cstddef>
#include <vector>

namespace dynsyslab::core {

struct KROMRolloutRequest {
    std::vector<std::complex<double>> eigenvalues;
    std::vector<std::complex<double>> modes;
    std::vector<std::complex<double>> initialCoordinates;
    std::vector<double> time;
    std::size_t stateDimension{};
};

struct KROMRollout {
    std::vector<std::complex<double>> lifted;
    std::vector<std::complex<double>> reconstructed;
    std::size_t samples{};
    std::size_t modes{};
    std::size_t stateDimension{};
};

[[nodiscard]] KROMRollout rolloutKROM(const KROMRolloutRequest& request);

}  // namespace dynsyslab::core
