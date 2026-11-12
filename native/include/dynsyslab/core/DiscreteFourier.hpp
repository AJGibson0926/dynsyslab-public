#pragma once
#include <complex>
#include <vector>
namespace dynsyslab::core {
// Unnormalized forward DFT; inverse includes 1/N. Arbitrary lengths, O(N log N).
[[nodiscard]] std::vector<std::complex<double>>
discreteFourierTransform(const std::vector<std::complex<double>> &input, bool inverse = false);
} // namespace dynsyslab::core
