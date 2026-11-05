#pragma once

#include <cstddef>
#include <vector>

namespace dynsyslab::numerics {

// Solves (X'X + ridge I) K = X'Y by complete orthogonal decomposition.
// All matrices are contiguous row-major runtime values.
[[nodiscard]] std::vector<double> solveNormalEquations(
    const std::vector<double>& matrix,std::size_t rows,std::size_t columns,
    const std::vector<double>& rightHandSide,std::size_t outputs,
    double ridge=0.0);

}  // namespace dynsyslab::numerics
