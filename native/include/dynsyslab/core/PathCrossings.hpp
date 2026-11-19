#pragma once
#include <cstddef>
#include <vector>

namespace dynsyslab::core {

// A transverse self-intersection of a sampled planar path: segment `first`
// (samples first..first+1) meets segment `second` at fractions u and v along
// them. first < second - 1, so neighbouring segments never count.
struct PathCrossing {
    std::size_t first{}, second{};
    double u{}, v{};
};

// All isolated, nonparallel self-intersections of the polyline (x[k], y[k]).
// Segments are half-open so a shared vertex is reported once; the final
// endpoint is included. Semantics match the miscellaneous-computes crossing
// kernel (mc_crossings_v1); the search uses a uniform grid, so long
// oscillating records cost about linear time. Results are ordered by
// (first, second).
[[nodiscard]] std::vector<PathCrossing> findPathCrossings(const std::vector<double>& x,
                                                          const std::vector<double>& y);

} // namespace dynsyslab::core
