#pragma once
// Observability of each eigenfunction of a K-ROM from a measurement: how well psi_n can
// be determined from the output over a window when every other eigenfunction is also
// unknown.
//
// With the finite-horizon modal gramian G = int_0^T e^{L* t} C* C e^{L t} dt (L = diag(lambda),
// one complex coordinate per eigenfunction) and white output noise of intensity V:
//   signal_n          = sqrt(G_nn)                 output RMS energy from a unit psi_n alone
//   distinguishable_n = 1 / sqrt((G^-1)_nn)        the part of psi_n's output that no
//                                                  combination of the other eigenfunctions
//                                                  can imitate (||P_perp O e_n||)
//   isolation_n       = distinguishable_n / signal_n  in [0, 1]
//   boundStd_n        = sqrt(V (G^-1)_nn)          Cramer-Rao bound on the std of psi_n-hat
//   pairBoundStd_n    = sqrt(V lambda_max(2x2 block of G^-1 over n and its conjugate partner)):
//                       the worst-direction bound for the pair the real output determines jointly
//   kalmanStd_n       = steady-state KLQE std of the conjugate-consistent estimate
//                       (psi_n + conj(psi_partner)) / 2 for process
//                       covariance W I on the split coordinates (optional; includes the dynamics)
// G^-1 uses a relative Tikhonov floor (regularization * max eig(G)), so numerically
// unobservable eigenfunctions get a large but finite bound rather than infinity.
#include <dynsyslab/core/KoopmanModel.hpp>

#include <complex>
#include <string>
#include <vector>

namespace dynsyslab::core {

struct EigenfunctionObservabilityOptions {
    double horizon{1.0};           // output window T (model time units)
    double noiseIntensity{1e-4};   // V
    double regularization{1e-12};  // relative Tikhonov floor on G
    double processCovariance{0.0}; // W > 0: also the KLQE steady-state std per eigenfunction
};

struct EigenfunctionObservability {
    std::size_t modes{};
    double horizon{}, noiseIntensity{};
    std::vector<std::complex<double>> eigenvalues;
    std::vector<std::size_t> partner;  // conjugate partner index (itself for a real mode)
    std::vector<double> signal, distinguishable, isolation, boundStd, pairBoundStd, kalmanStd;
    std::vector<std::complex<double>> gramian;  // modes x modes, Hermitian
};

[[nodiscard]] EigenfunctionObservability eigenfunctionObservability(const KoopmanModel& model,
    const std::string& measurement, const EigenfunctionObservabilityOptions& options = {});

}  // namespace dynsyslab::core
