#pragma once
// State samples conditioned on a path, for fitting an approximate input map where a controller
// will use it:
//   resampleByArcLength: `count` states equally spaced in arc length along a piecewise-linear
//     path (uniform in distance, not time, so slow segments are not over-weighted);
//   tubeStates: each path state plus `perPoint` states displaced by `radius` times the local
//     axis scale in random directions (seeded), a tube around the path.
#include <cstddef>
#include <cstdint>
#include <vector>

namespace dynsyslab::core {

[[nodiscard]] std::vector<double> resampleByArcLength(const std::vector<double>& states, std::size_t samples,
                                                      std::size_t dimension, std::size_t count);
[[nodiscard]] std::vector<double> tubeStates(const std::vector<double>& path, std::size_t samples,
                                             std::size_t dimension, double radius, std::size_t perPoint,
                                             std::uint64_t seed);

}  // namespace dynsyslab::core
