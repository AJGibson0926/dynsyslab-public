#pragma once
// Internal pieces shared by the Koopman SINDy fits (SINDy.cpp, SindyEigenvalues.cpp).
#include <complex>
#include <cstddef>
#include <optional>
#include <vector>

namespace dynsyslab::core::sindy_detail {

// Per-axis mean and standard deviation (0 and 1 when `normalize` is false).
void normalization(const std::vector<double>& states, std::size_t samples, std::size_t dimension, bool normalize,
                   std::vector<double>& center, std::vector<double>& scale);
// Theta: the polynomial library at each sample (samples x features), on normalized variables.
std::vector<double> featureMatrix(const std::vector<double>& states, std::size_t samples, std::size_t dimension,
                                  const std::vector<std::size_t>& exponents, const std::vector<double>& center,
                                  const std::vector<double>& scale);
// Second-order, nonuniform-step finite differences along the rows.
std::vector<double> finiteDifference(const std::vector<double>& values, std::size_t samples, std::size_t columns,
                                     const std::vector<double>& time);
// Gamma: the rate of every library function along the dynamics, grad(theta)·xdot (samples x features).
std::vector<double> featureDerivative(const std::vector<double>& states, const std::vector<double>& velocities,
                                      std::size_t samples, std::size_t dimension,
                                      const std::vector<std::size_t>& exponents, const std::vector<double>& center,
                                      const std::vector<double>& scale);
// Sequentially thresholded smallest right singular vector of a tall complex matrix (unit norm).
std::vector<std::complex<double>> sparseNullVector(const std::vector<std::complex<double>>& source, std::size_t rows,
                                                   std::size_t columns, double threshold,
                                                   std::size_t maximumIterations);

// Real feature Gram blocks, formed once per fit: A = ThetaᵀTheta, B = GammaᵀTheta,
// C = GammaᵀGamma (all features x features, row-major).
struct GramBlocks { std::vector<double> A, B, C; std::size_t n{}; };
GramBlocks gramBlocks(const std::vector<double>& theta, const std::vector<double>& gamma, std::size_t rows,
                      std::size_t n);
// The thresholded null vector of G(lambda) = (Gamma - lambda Theta)^H (Gamma - lambda Theta); nullopt
// when the Gram system cannot be trusted (callers fall back to sparseNullVector).
std::optional<std::vector<std::complex<double>>> sparseNullVectorGram(const GramBlocks& g, std::complex<double> lambda,
                                                                      double threshold,
                                                                      std::size_t maximumIterations);
bool conjugatePartner(std::complex<double> left, std::complex<double> right);
// Dense complex solve (Gaussian elimination with partial pivoting).
std::vector<std::complex<double>> solveDense(std::vector<std::complex<double>> matrix,
                                             std::vector<std::complex<double>> rhs, std::size_t n);
// The anchored fit's coefficients: the thresholded solution of
// (Theta_aᵀTheta_a + w G(lambda)) xi = Theta_aᵀ psi, w = weight tr(Theta_aᵀTheta_a) / tr(G(lambda)).
std::vector<std::complex<double>> anchoredVector(const std::vector<double>& anchorGram, const GramBlocks& g,
                                                 std::complex<double> lambda,
                                                 const std::vector<std::complex<double>>& rhs, double weight,
                                                 double threshold, std::size_t maximumIterations);

}  // namespace dynsyslab::core::sindy_detail
