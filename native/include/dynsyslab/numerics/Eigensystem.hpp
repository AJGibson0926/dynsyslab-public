#pragma once

#include <complex>
#include <cstddef>
#include <vector>

namespace dynsyslab::numerics {

struct RealEigensystem {
    std::vector<std::complex<double>> eigenvalues;
    std::vector<std::complex<double>> rightEigenvectors; // row-major
    std::size_t dimension{};
};

// A = U T U*, with U unitary and T upper triangular; arrays are row-major.
// Retains the full triangular factor, including nonnormal/repeated-root terms.
struct ComplexSchurDecomposition {
    std::vector<std::complex<double>> unitary,triangular;
    std::size_t dimension{};
};
[[nodiscard]] ComplexSchurDecomposition complexSchur(
    const std::vector<double>& matrix,std::size_t dimension);

[[nodiscard]] RealEigensystem eigensystem(
    const std::vector<double>& matrix,std::size_t dimension,
    bool computeEigenvectors=true);

}  // namespace dynsyslab::numerics
