#pragma once
#include <cstddef>
#include <vector>

namespace dynsyslab::core {

// Error-bounded vertex selection for drawing or storing a long sampled path
// (after the trajectory-resolution side-car, tools/experiments/
// trajectory_resolution/vertices.cpp, generalized to any number of columns).
// The first and last samples and each column's global extrema are kept; every
// dropped sample lies within `fraction` of each column's range from the
// straight segment between the kept samples around it, at the same fraction of
// the index interval. Failing segments are split at their midpoint, so fast
// excursions keep their vertices and smooth stretches are thinned.
// values: rows x width, row-major. Returns increasing row indices.
[[nodiscard]] std::vector<std::size_t> boundedPathVertices(const std::vector<double>& values,
    std::size_t rows, std::size_t width, double fraction);

// The same selection with the tightest tolerance (starting at `fraction`,
// doubling) that keeps at most `maximum` vertices. Falls back to an even
// stride when the values are not finite.
[[nodiscard]] std::vector<std::size_t> simplifiedPathIndices(const std::vector<double>& values,
    std::size_t rows, std::size_t width, std::size_t maximum, double fraction = 1e-4);

} // namespace dynsyslab::core
