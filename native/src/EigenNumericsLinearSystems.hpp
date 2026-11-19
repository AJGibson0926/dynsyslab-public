#pragma once
// Private Eigen helpers shared by the linear-systems toolkit sources.
#include <dynsyslab/core/StateSpace.hpp>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <complex>
#include <vector>

namespace dynsyslab::core::linear_detail {

using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using CMatrix = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>;
using Vector = Eigen::VectorXd;

[[nodiscard]] Matrix mapped(const std::vector<double>& values, std::size_t rows,
                            std::size_t columns, const char* name);
[[nodiscard]] std::vector<double> flat(const Matrix& matrix);

struct Realization { Matrix A, B, C, D; };
[[nodiscard]] Realization realization(const StateSpaceModel& model);
[[nodiscard]] StateSpaceModel fromRealization(const Realization& r, TimeDomain domain,
                                              double sampleTime);

[[nodiscard]] std::vector<std::complex<double>> eigenvalues(const Matrix& matrix);
[[nodiscard]] bool stable(const Matrix& A, bool discrete, double margin = 0.0);

// The one Lyapunov/Stein solver in DynSysLab (Bartels-Stewart on the complex
// Schur form): A X + X A^T + Q = 0 (continuous) or A X A^T - X + Q = 0 (discrete).
// numerics' infinite observability gramians call it with A^T.
[[nodiscard]] Matrix lyapunov(const Matrix& A, const Matrix& Q, bool discrete);
// The one finite-horizon continuous gramian: W(T) = int_0^T e^{A t} Q e^{A^T t} dt by
// Van Loan's block exponential on a short step, then horizon doubling
// W(2h) = W(h) + e^{A h} W(h) e^{A^T h} (no overflow for stiff stable A).
[[nodiscard]] Matrix finiteContinuousGramian(const Matrix& A, const Matrix& Q, double horizon);
// Square-root factor L (n x n) with L L^T = W for symmetric PSD W; negative
// rounding eigenvalues are clamped to zero.
[[nodiscard]] Matrix psdFactor(const Matrix& W);
// Orthonormal basis for the column space with relative singular-value tolerance.
[[nodiscard]] Matrix orthonormalBasis(const Matrix& M, double relativeTolerance);
// Orthonormal basis of the orthogonal complement of range(Q) (Q orthonormal).
[[nodiscard]] Matrix orthogonalComplement(const Matrix& Q, Eigen::Index dimension);

}  // namespace dynsyslab::core::linear_detail
