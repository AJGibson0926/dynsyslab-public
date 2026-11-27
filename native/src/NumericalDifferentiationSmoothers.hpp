#pragma once
// Model-based derivative estimators of NumericalDifferentiation.hpp, split from its dispatcher:
// the Rauch-Tung-Striebel smoother of an integrated white-noise model (Kalman 1960; Rauch, Tung &
// Striebel 1965), direct 1-D total-variation denoising (Condat 2013), and the integral-form
// regularized derivative (Chartrand 2011).
#include <dynsyslab/core/NumericalDifferentiation.hpp>

#include <vector>

namespace dynsyslab::core::differentiation_detail {
double rmsAboutMean(const std::vector<double>& y);
std::vector<double> centralDifference(const std::vector<double>& y, double dt);
// The integrated white-noise model of order d (state y, y', ..., y^(d-1)): transition F, noise Q.
void integratedModel(std::size_t d, double dt, double q, std::vector<double>& F, std::vector<double>& Q);
double kalmanQ(const DerivativeSpec& spec, double r, double dt, std::size_t d);
std::vector<double> matmul(const std::vector<double>& A, const std::vector<double>& B, std::size_t d);
std::vector<double> transpose(const std::vector<double>& A, std::size_t d);
DerivativeResult kalmanSmoother(const std::vector<double>& y, double dt, const DerivativeSpec& spec);
std::vector<double> tvDenoise(const std::vector<double>& input, double lambda);
DerivativeResult integralRegularized(const std::vector<double>& y, double dt, const DerivativeSpec& spec, bool tv);
}  // namespace dynsyslab::core::differentiation_detail
