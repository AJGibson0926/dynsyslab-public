#pragma once
// Candidate libraries for vector-field and map regression beyond the polynomial and mixed
// families of DynamicsLibrary.hpp (docs/architecture/TRAJECTORY_KERNELS.md, "Model libraries").
// Every family is frozen on training coordinates and returned as analytic expressions, so fits,
// derivatives and rollouts use the shared expression engine. Ported from the miscellaneous-computes
// side-car (features.prepare_library / library). With z = x / sqrt(1 + (x/3)^2) (bounded coordinates):
//
//   polynomial    1 and the monomials of x up to `degree`
//   bounded       1, x and the monomials of z
//   rational      bounded plus w^2/(w^2 + x_j^2) and w x_j/(w^2 + x_j^2) per width w (regularized
//                 reciprocals; rational SINDy terms, Kaheman, Kutz and Brunton 2020, Proc. R. Soc. A
//                 476:20200279)
//   localized     bounded plus a normalized partition of unity over `centers` training points
//                 (farthest-point sampling, Gonzalez 1985, Theor. Comput. Sci. 38:293-306), weights
//                 (1 + |z - c|^2 / s_c)^-2 / sum, times [1, z] (Shepard 1968 normalization)
//   gated         g(E) x_j with E the mean square of x and g in {1, E/(E + s), s = .04, .4, 4, 16}
//   paired_gated  gated plus the per-pair gates p/(1 + p) of consecutive coordinate pairs
//   harmonic      x and, for k = 1..harmonics, a(x) Re u^k and a(x) Im u^k with
//                 u = (x_1 + i x_2)/sqrt(|.|^2 + .015) and amplitude gates a (a phase-amplitude
//                 basis for a dominant modal pair)
// Candidate regressions remain sequentially thresholded least squares (Brunton, Proctor and Kutz
// 2016, PNAS 113:3932-3937).
#include <dynsyslab/core/Observable.hpp>
#include <dynsyslab/core/TimeSeries.hpp>

#include <cstddef>
#include <string>
#include <vector>

namespace dynsyslab::core {

struct ModelLibraryOptions {
    std::size_t degree{2};
    std::size_t centers{24};
    std::size_t harmonics{5};
    std::vector<double> rationalWidths{0.25, 1.0, 4.0};
};

[[nodiscard]] std::vector<std::string> modelLibraryFamilies();

// The library of `family` frozen on `training` (samples x coordinates).
[[nodiscard]] Observable modelLibrary(const std::string& family, const TimeSeries& training, const ModelLibraryOptions& options = {});

}  // namespace dynsyslab::core
