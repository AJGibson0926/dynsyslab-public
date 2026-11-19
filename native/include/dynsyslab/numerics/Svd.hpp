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

// Complex SVD M = U diag(s) V^H (thin): singular values descending; right vectors
// columns x min(rows, columns), row-major, singular vector k in column k.
struct ComplexSvd {
    std::vector<double> singularValues;
    std::vector<std::complex<double>> right;
    std::size_t rows{}, columns{}, rank{};
};
[[nodiscard]] ComplexSvd complexSingularValueDecomposition(
    const std::vector<std::complex<double>>& matrix,std::size_t rows,std::size_t columns);

[[nodiscard]] std::vector<double> singularValues(
    const std::vector<std::complex<double>>& matrix,std::size_t rows,
    std::size_t columns);

}  // namespace dynsyslab::numerics
