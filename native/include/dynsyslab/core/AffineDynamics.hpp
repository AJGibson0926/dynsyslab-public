#pragma once
// Affine vector fields x' = A x + c: their equations from a fitted model and their exact flow
// (docs/architecture/TRAJECTORY_KERNELS.md).
//
//   affineEquations     A = Df(0), c = f(0) of a vector field, verified affine: f must equal
//                       A x + c at probe points across `scale` (otherwise the model has active
//                       nonlinear terms and no A, c form).
//   affinePropagation   x(t) from the augmented generator [[A, c], [0, 0]]: exp(t M) [x0; 1]
//                       (Van Loan 1978, IEEE Trans. Autom. Control 23:395-404), on a uniform clock
//                       by repeated doubling of one step's exponential.
#include <dynsyslab/core/Observable.hpp>

#include <cstddef>
#include <vector>

namespace dynsyslab::core {

struct AffineEquations {
    std::size_t dimension{};
    std::vector<double> A;  // dimension x dimension, row-major
    std::vector<double> c;  // dimension
};

[[nodiscard]] AffineEquations affineEquations(const Observable& field, double scale = 1.0, double tolerance = 1e-9);

// States (times x dimension) on `times` (increasing, times[0] is the initial time).
[[nodiscard]] std::vector<double> affinePropagation(const AffineEquations& equations, const std::vector<double>& initial,
                                                    const std::vector<double>& times);

}  // namespace dynsyslab::core
