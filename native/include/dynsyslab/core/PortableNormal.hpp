#pragma once

#include <cstdint>
#include <random>

namespace dynsyslab::core {

// Identifier of the random stream below. Archives record it; change it if the stream changes.
inline constexpr const char* kPortableNormalStream = "mt19937_64+as241-inverse-cdf/v1";

// Standard normal quantile by Wichura's PPND16 (Algorithm AS 241, Applied Statistics 37(3),
// 477-484, 1988). Uses only IEEE +, -, *, / and sqrt, with a fixed-series logarithm, so the
// value does not depend on a compiler's math library. p must lie in (0, 1).
[[nodiscard]] double portableNormalQuantile(double p);

// One standard normal from one 64-bit Mersenne Twister word (std::mt19937_64 is fully specified
// by the C++ standard, unlike std::normal_distribution). The top 53 bits select p in (0, 1).
[[nodiscard]] double portableNormal(std::mt19937_64& generator);

}  // namespace dynsyslab::core
