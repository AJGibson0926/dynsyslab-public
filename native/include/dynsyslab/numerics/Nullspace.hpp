#pragma once

#include <complex>
#include <cstddef>
#include <vector>

namespace dynsyslab::numerics {

[[nodiscard]] std::vector<std::complex<double>> smallestRightSingularVector(
    const std::vector<std::complex<double>>& matrix,std::size_t rows,
    std::size_t columns);

}  // namespace dynsyslab::numerics
