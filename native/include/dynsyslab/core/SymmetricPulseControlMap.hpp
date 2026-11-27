#pragma once
// B_psi(x) estimated from paired opposite-command first responses, without the analytic map.
// A pair starts both runs at one anchor x_0 and holds -u and +u for one sample dt; since
// psi(x_+) - psi(x_-) = 2 u dt B_psi(x_0) + O(dt^2) (the drift cancels to first order),
//   B_psi(x_0) ~= (psi(x_+(dt)) - psi(x_-(dt))) / ((u_+ - u_-) dt),
// a central difference in the input (Korda & Mezic 2018, "Linear predictors for nonlinear
// dynamical systems: Koopman operator meets model predictive control", Automatica 93:149-160,
// for the input-affine lifted model it estimates). The samples then fit a constant map (their
// mean) or a bilinear map b_0 + sum_i b_i psi_i(x_0) by least squares.
#include <dynsyslab/core/KoopmanModel.hpp>

#include <complex>
#include <cstddef>
#include <string>
#include <vector>

namespace dynsyslab::core {

struct SymmetricPulseRun {
    std::size_t pairId{};
    int polarity{};                // -1 or +1
    std::vector<double> anchor;    // x_0 (model coordinates)
    std::vector<double> response;  // x(dt)
    double input{};                // the held command
    double sampleTime{};
};

struct SymmetricPulseMapSamples {
    std::vector<double> anchors;                // pairs x state
    std::vector<std::complex<double>> maps;     // pairs x modes (one input)
    std::vector<double> responseNorms;          // |psi_+ - psi_-| per pair
    std::size_t pairs{}, modes{}, stateDimension{};
};

[[nodiscard]] SymmetricPulseMapSamples estimateSymmetricPulseControlMap(
    const KoopmanModel& model, const std::vector<SymmetricPulseRun>& runs);

struct SampledMapFit {
    ApproximateControlMap map;
    double relativeResidual{};
    double complexRmse{};
};

[[nodiscard]] SampledMapFit fitConstantControlMapSamples(const KoopmanModel& model,
                                                         const SymmetricPulseMapSamples& samples);
[[nodiscard]] SampledMapFit fitBilinearControlMapSamples(const KoopmanModel& model,
                                                         const SymmetricPulseMapSamples& samples, double ridge = 0.0);

}  // namespace dynsyslab::core
