#pragma once
// Post-fit editing of expansion eigenfunctions psi_j(q) = sum_k c_jk theta_k(q) of a K-ROM:
//
//  refitAtMultiplier: the coefficients that best satisfy psi(q_{k+1}) = mu psi(q_k) along sampled
//    runs for a fixed multiplier mu (the eigenvalue is kept), as the right singular vector of
//    D = Theta(q_{k+1}) - mu Theta(q_k) with the smallest singular value (the null-space
//    eigenfunction fit of SINDy-Koopman: Kaiser, Kutz & Brunton 2021, "Data-driven discovery of
//    Koopman eigenfunctions for control", Mach. Learn.: Sci. Technol. 2:035023).
//  gaugeFactor: the complex scale alpha minimising ||alpha psi - psi_ref|| on samples, so a refit
//    or re-identified eigenfunction keeps the reference's normalisation and phase (an
//    eigenfunction is defined up to a constant factor).
//  pruneCoefficients: zero terms with |c_k| < threshold * max_k |c_k|.
#include <complex>
#include <cstddef>
#include <vector>

namespace dynsyslab::core {

// thetaNow, thetaNext: pairs x terms (row-major dictionary values). Returns unit-norm coefficients.
[[nodiscard]] std::vector<std::complex<double>> refitAtMultiplier(
    const std::vector<double>& thetaNow, const std::vector<double>& thetaNext,
    std::size_t pairs, std::size_t terms, std::complex<double> multiplier, double* residual = nullptr);

[[nodiscard]] std::complex<double> gaugeFactor(const std::vector<std::complex<double>>& values,
                                               const std::vector<std::complex<double>>& reference);

// Returns the number of terms set to zero.
std::size_t pruneCoefficients(std::vector<std::complex<double>>& coefficients, double threshold);

}  // namespace dynsyslab::core
