#pragma once

#include <complex>
#include <cstddef>
#include <vector>

namespace dynsyslab::numerics {

enum class LeastSquaresMethod {
    columnPivotedQr,
    divideAndConquerSvd
};

struct LeastSquaresOptions {
    LeastSquaresMethod method{LeastSquaresMethod::columnPivotedQr};
    double ridge{};
    double singularValueTolerance{};
};

[[nodiscard]] std::vector<double> solveLeastSquares(
    const std::vector<double>& matrix,std::size_t rows,std::size_t columns,
    const std::vector<double>& rightHandSide,std::size_t outputs,
    const LeastSquaresOptions& options={});

[[nodiscard]] std::vector<std::complex<double>> solveLeastSquares(
    const std::vector<std::complex<double>>& matrix,std::size_t rows,
    std::size_t columns,
    const std::vector<std::complex<double>>& rightHandSide,
    std::size_t outputs,const LeastSquaresOptions& options={});

}  // namespace dynsyslab::numerics
