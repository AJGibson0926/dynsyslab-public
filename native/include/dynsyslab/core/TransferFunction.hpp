#pragma once
// Transfer functions of one input/output channel of an LTI realization (e.g. a
// K-ROM's split realization from koopmanStateSpace): pole/zero/gain form,
// polynomial form, and Bode data with unwrapped phase.
//
// Everything is on the realization as given: poles are eig(A); zeros are the
// invariant zeros of the channel (C_i, A, B_j, D_ij). Uncontrollable or
// unobservable modes therefore appear as near-coincident pole/zero pairs
// (cancellations) rather than being removed.
#include <dynsyslab/core/StateSpace.hpp>

#include <complex>
#include <cstddef>
#include <vector>

namespace dynsyslab::core {

struct TransferFunction {
    std::size_t output{}, input{};
    std::vector<std::complex<double>> poles, zeros;
    double gain{};                    // G(s) = gain * prod(s - z) / prod(s - p)
    std::vector<double> numerator;    // descending powers, length zeros + 1 (gain * monic zero polynomial)
    std::vector<double> denominator;  // descending powers, monic, length poles + 1
    int relativeDegree{};             // poles - zeros
    double directGain{};              // D_ij
    std::vector<std::pair<std::size_t, std::size_t>> cancellations;  // (pole, zero) pairs closer than the tolerance
};

// The single-input single-output channel (output, input) of the model.
[[nodiscard]] StateSpaceModel channel(const StateSpaceModel& model, std::size_t output, std::size_t input);

[[nodiscard]] TransferFunction transferFunction(const StateSpaceModel& model, std::size_t output, std::size_t input,
                                                double cancellationTolerance = 1e-6);

// Monic polynomial with the given roots (descending powers, real parts; the
// roots must come in conjugate pairs).
[[nodiscard]] std::vector<double> polynomialFromRoots(const std::vector<std::complex<double>>& roots);

struct BodeData {
    std::vector<double> frequencies;  // rad per model time unit
    std::vector<double> magnitudeDb;  // 20 log10 |G|
    std::vector<double> phaseDeg;     // unwrapped along the grid, starting in (-180, 180]
};
[[nodiscard]] BodeData bode(const StateSpaceModel& model, std::size_t output, std::size_t input,
                            const std::vector<double>& frequencies);
// Unwrapped phase in degrees of a complex response sampled along a frequency grid.
[[nodiscard]] std::vector<double> unwrappedPhaseDeg(const std::vector<std::complex<double>>& response);

}  // namespace dynsyslab::core
