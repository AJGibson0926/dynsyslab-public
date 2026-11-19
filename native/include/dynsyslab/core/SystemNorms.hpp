#pragma once
// Input/output characteristics of an LTI realization: frequency response,
// H2 and H-infinity norms, transmission zeros, DC gain, modal damping and
// SISO stability margins.
#include <dynsyslab/core/StateSpace.hpp>

#include <complex>
#include <cstddef>
#include <vector>

namespace dynsyslab::core {

struct FrequencyResponse {
    std::vector<double> frequencies;                 // rad per model time unit
    std::vector<std::complex<double>> response;      // frequency x outputs x inputs
    std::vector<double> singularValues;              // frequency x min(outputs, inputs), descending
    std::size_t outputs{}, inputs{};
};
// Continuous: G(j w). Discrete: G(exp(j w dt)); frequencies above pi/dt alias.
[[nodiscard]] FrequencyResponse frequencyResponse(const StateSpaceModel& model,
    const std::vector<double>& frequencies);
// Logarithmic grid spanning two decades beyond the pole/zero magnitudes
// (capped at the Nyquist frequency for discrete models).
[[nodiscard]] std::vector<double> logarithmicFrequencyGrid(const StateSpaceModel& model,
    std::size_t count = 400);

// sqrt(trace(C Wc C^T)) (+ trace(D D^T) for discrete). Infinite for unstable
// models or a continuous model with nonzero D.
[[nodiscard]] double h2Norm(const StateSpaceModel& model);

struct HInfinityNorm {
    double value{};          // sup over frequency of the largest singular value
    double peakFrequency{};  // rad per model time unit
    std::size_t iterations{};
};
// Boyd-Balakrishnan / Bruinsma-Steinbuch Hamiltonian level-set iteration.
// Discrete models are mapped through the bilinear transform. Infinite for
// unstable models.
[[nodiscard]] HInfinityNorm hInfinityNorm(const StateSpaceModel& model, double relativeTolerance = 1.0e-9);

// Invariant zeros of a square system: finite generalized eigenvalues of the
// Rosenbrock pencil [[A, B], [C, D]] - s [[I, 0], [0, 0]]. Non-square systems
// are rejected explicitly.
[[nodiscard]] std::vector<std::complex<double>> transmissionZeros(const StateSpaceModel& model);

// G(0) (continuous) or G(1) (discrete), outputs x inputs.
[[nodiscard]] std::vector<double> dcGain(const StateSpaceModel& model);

struct ModalDamping {
    std::complex<double> pole;          // continuous-equivalent pole for discrete models
    double naturalFrequency{}, dampingRatio{}, timeConstant{};
};
[[nodiscard]] std::vector<ModalDamping> modalDamping(const StateSpaceModel& model);

// Classical margins of the negative-feedback loop around a SISO open loop L.
struct StabilityMargins {
    double gainMargin{};        // linear factor; infinite if the phase never crosses -180 deg
    double phaseMargin{};       // degrees; infinite if |L| never crosses 1
    double phaseCrossover{};    // rad per time unit
    double gainCrossover{};
};
[[nodiscard]] StabilityMargins stabilityMargins(const StateSpaceModel& openLoop, std::size_t gridPoints = 4000);

}  // namespace dynsyslab::core
