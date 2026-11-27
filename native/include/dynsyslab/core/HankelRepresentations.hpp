#pragma once
// Hankel DMD's eigenfunctions as a function of series windows (the "time-delay" representation):
// psi(h) = W^-1 U^T h = (Phi^H Phi)^-1 Phi^H h, with Phi = U W the Hankel modes
// (Arbabi & Mezic, SIAM J. Appl. Dyn. Syst. 16, 2017; Brunton et al., Nat. Commun. 8, 2017).
// It is an ordinary FunctionExpansion: one coordinate term per window row, on the window domain
// (FunctionDomain::seriesWindow), so an Observable evaluates it along any series by cutting the
// window at each sample. Row j of T = G^-1 Phi^H (G = Phi^H Phi) is eigenfunction j's own linear map.
#include <dynsyslab/core/FunctionExpansion.hpp>
#include <dynsyslab/core/Hankel.hpp>
#include <dynsyslab/core/Json.hpp>

#include <complex>
#include <memory>
#include <vector>

namespace dynsyslab::core {

// From the Hankel modes Phi (rows x rank) and their Gram matrix G (rank x rank): all `rank` rows of T.
[[nodiscard]] std::shared_ptr<const FunctionExpansion> windowMapFromModes(const HankelWindowLayout& layout,
    const std::vector<std::complex<double>>& modes, const std::vector<std::complex<double>>& gram, std::size_t rank);
// From the map T itself (modes x rows, applied to raw windows).
[[nodiscard]] std::shared_ptr<const FunctionExpansion> windowMapFromMap(const HankelWindowLayout& layout,
    const std::vector<std::complex<double>>& map, std::size_t modes);
// T (modes x rows) of a window map, for archives; the layout is the expansion's domain window.
[[nodiscard]] std::vector<std::complex<double>> windowMapMatrix(const FunctionExpansion& windowMap);
// The archive description of a window map: its layout and form.
[[nodiscard]] json::Value windowMapDescription(const HankelWindowLayout& layout);

}  // namespace dynsyslab::core
