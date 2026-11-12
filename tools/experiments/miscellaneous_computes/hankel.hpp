#pragma once
#include <cstddef>

// Add H H^T for one trajectory. H columns are complete causal windows.
// Delay-major, channel-minor order. Both input and Gram array are row-major.
// The caller zeros Gram once and may accumulate additional independent records.
// Returns 0 for success, 1 for invalid dimensions/data, 2 for an allocation error.
extern "C" int mc_hankel_gram_v1(const double* observations, std::size_t samples,
    std::size_t channels, std::size_t taps, double* gram) noexcept;
