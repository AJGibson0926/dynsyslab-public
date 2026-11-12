#pragma once

#include <complex>
#include <cstddef>
#include <vector>

namespace dynsyslab::numerics {

struct RealSvd {
    std::vector<double> left;             // rows x rank
    std::vector<double> singularValues;   // rank
    std::vector<double> rightTranspose;   // rank x columns
    std::size_t rows{};
    std::size_t columns{};
    std::size_t rank{};
};

[[nodiscard]] RealSvd singularValueDecomposition(
    const std::vector<double>& matrix,std::size_t rows,std::size_t columns,
    std::size_t requestedRank=0,double relativeTolerance=0.0);

[[nodiscard]] std::vector<double> singularValues(
    const std::vector<std::complex<double>>& matrix,std::size_t rows,
    std::size_t columns);

}  // namespace dynsyslab::numerics
