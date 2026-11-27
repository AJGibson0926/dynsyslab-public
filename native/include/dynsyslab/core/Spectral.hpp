#pragma once

#include <complex>
#include <cstddef>
#include <string>
#include <vector>

#include <dynsyslab/core/Observable.hpp>
#include <dynsyslab/core/Spectrum.hpp>
#include <dynsyslab/core/Hankel.hpp>

namespace dynsyslab::core {

struct SVDResult {
    std::vector<double> left;
    std::vector<double> singularValues;
    std::vector<double> right;
    std::size_t rows{};
    std::size_t columns{};
    std::size_t rank{};
};

struct DMDOptions {
    std::size_t rank{};              // zero selects numerical rank
    double singularValueTolerance{}; // zero selects the numerical threshold
    bool exactModes{true};
    // Gavish & Donoho (2014) optimal hard threshold for unknown white noise:
    // keep singular values above omega(beta) * median(sigma), beta = min(m,n)/max(m,n).
    // Applied before `rank` (which still caps the result).
    bool optimalHardThreshold{false};
};

struct DMDResult {
    Spectrum spectrum;  // the sampled multipliers mu and their principal rates log(mu) / deltaTime
    std::vector<std::complex<double>> modes;       // featureDimension x rank
    std::vector<std::complex<double>> coordinates; // samples x rank
    std::vector<std::complex<double>> operatorMatrix;
    std::vector<std::complex<double>> eigenvectors; // rank x rank W (projected DMD with coordinates)
    std::size_t featureDimension{};
    std::size_t samples{};
    std::size_t rank{};
    double deltaTime{};
};

struct HankelDMDOptions : DMDOptions {
    std::size_t matrixSize{2000};
    std::size_t samplingWindow{20000};
    std::size_t sampleStride{1};
    std::size_t observableIndex{};
    bool terminalAmplitudeEnergySort{true};
    // Forward windows retain the historical start alignment by default.
    // End alignment permits causal evaluation once a full window is available.
    bool alignToWindowEnd{false};
    std::size_t snapshotStride{1}; // spacing between training window starts; paired shift remains one sample
    std::vector<double> channelScales; // optional positive divisors; no centering
};

struct HankelDMDResult : DMDResult {
    std::vector<std::complex<double>> modalAmplitudes;
    std::vector<std::size_t> sourceModeIndices;
    std::size_t matrixRows{};
    std::size_t matrixColumns{};
    std::size_t sourceSamples{};
    std::string observableName;
    std::vector<std::string> observableNames;
    std::vector<std::complex<double>> reconstructionMatrix; // channels x rank
    std::vector<double> singularValues;
    std::size_t delayCount{};
    std::size_t referenceDelay{};
    std::vector<std::size_t> coordinateSampleOffsets; // original source sample indices
    // Phi^H Phi = W^H W (rank x rank, U orthonormal). With the modes it gives the exact
    // window-to-eigenfunction map psi(h) = W^-1 U^T h = G^-1 Phi^H h.
    std::vector<std::complex<double>> modeGram;
};

[[nodiscard]] SVDResult singularValueDecomposition(
    const std::vector<double>& rowMajorValues, std::size_t rows,
    std::size_t columns, std::size_t rank = 0);

[[nodiscard]] DMDResult dynamicModeDecomposition(
    const std::vector<double>& rowMajorStates, std::size_t samples,
    std::size_t stateDimension, double deltaTime,
    const DMDOptions& options = {});

[[nodiscard]] HankelDMDResult hankelDynamicModeDecomposition(
    const std::vector<double>& rowMajorStates, std::size_t samples,
    std::size_t stateDimension, double deltaTime,
    const HankelDMDOptions& options = {});

// Fit an already constructed pair. Names describe the expanded input channels.
[[nodiscard]] HankelDMDResult fitHankelDMD(HankelMatrix matrix,
    const std::vector<std::string>& names, double sourceDeltaTime, const HankelDMDOptions& options = {});

// One shared decomposition, with delay-major rows: delay * channels + channel.
// Values are samples x channels. Every channel receives its own reconstruction
// row in observableNames order, in the original (unscaled) measurement units.
[[nodiscard]] HankelDMDResult hankelDynamicModeDecomposition(
    const std::vector<double>& rowMajorMeasurements, std::size_t samples,
    const std::vector<std::string>& observableNames, double deltaTime,
    const HankelDMDOptions& options = {});

[[nodiscard]] HankelDMDResult hankelDynamicModeDecomposition(
    const std::vector<double>& rowMajorStates, std::size_t samples,
    std::size_t stateDimension, double deltaTime,
    const std::vector<Observable>& observables,
    const HankelDMDOptions& options = {});

// H-DMD is defined on a scalar measurement, not intrinsically on a state
// component. This overload evaluates any real scalar observable first.
[[nodiscard]] HankelDMDResult hankelDynamicModeDecomposition(
    const std::vector<double>& rowMajorStates, std::size_t samples,
    std::size_t stateDimension, double deltaTime, const Observable& observable,
    const HankelDMDOptions& options = {});

}  // namespace dynsyslab::core
