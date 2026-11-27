#pragma once
// Distance between two finite eigenvalue sets (H-DMD spectra of different observables).
#include <complex>
#include <vector>

namespace dynsyslab::core {

// Symmetric nearest-neighbour RMS distance:
//   sqrt( (mean_i min_j |a_i - b_j|^2 + mean_j min_i |a_i - b_j|^2) / 2 ),
// the mean of the two directed (chamfer) distances, zero exactly when the sets coincide. NaN when
// either set is empty. (Borgefors, IEEE Trans. PAMI 10(6), 849, 1988, for chamfer matching.)
[[nodiscard]] double symmetricNearestRms(const std::vector<std::complex<double>>& a,
                                         const std::vector<std::complex<double>>& b);

}  // namespace dynsyslab::core
