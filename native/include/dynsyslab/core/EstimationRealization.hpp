#pragma once
#include <dynsyslab/core/KoopmanModel.hpp>
namespace dynsyslab::core {
struct RealModalRealization {
    KoopmanModel model;
    std::vector<std::size_t> sourceModes, conjugateModes; // self for real modes
    std::vector<std::complex<double>> divisors; // normalized psi = source psi / divisor
};
// Convert selected real modes/conjugate pairs from a sampled H-DMD KROM to
// independent RMS-normalized coordinates. Keeps the exact retained real-output
// reconstruction, including phase conventions; no output coefficients are fit.
// Positive-imaginary representatives are conventional, but either sign works.
// Pair relationships are verified on the archived training coordinates.
// A nominally real rate may be snapped to real within realEigenvalueTolerance
// (relative to max(1,abs(rate))); its sampled coordinate must also pass the
// real-up-to-phase check. Zero tolerance requires an exactly real rate.
[[nodiscard]] RealModalRealization makeRealModalRealization(
    const KoopmanModel& parent,const std::vector<std::size_t>& representatives,
    double relativeTolerance=1e-6,double realEigenvalueTolerance=1e-10);
} // namespace dynsyslab::core
