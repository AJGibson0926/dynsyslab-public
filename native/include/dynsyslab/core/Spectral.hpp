#pragma once

#include <complex>
#include <cstddef>
#include <string>
#include <vector>

#include <dynsyslab/core/Observable.hpp>

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
};

struct DMDResult {
    std::vector<std::complex<double>> discreteEigenvalues;
    std::vector<std::complex<double>> continuousEigenvalues;
    std::vector<std::complex<double>> modes;       // featureDimension x rank
    std::vector<std::complex<double>> coordinates; // samples x rank
    std::vector<std::complex<double>> operatorMatrix;
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
};

struct HankelDMDResult : DMDResult {
    std::vector<std::complex<double>> modalAmplitudes;
    std::vector<std::size_t> sourceModeIndices;
    std::size_t matrixRows{};
    std::size_t matrixColumns{};
    std::size_t sourceSamples{};
    std::string observableName;
    std::vector<std::complex<double>> reconstructionMatrix; // 1 x rank Koopman modes
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

// H-DMD is defined on a scalar measurement, not intrinsically on a state
// component. This overload evaluates any real scalar observable first.
[[nodiscard]] HankelDMDResult hankelDynamicModeDecomposition(
    const std::vector<double>& rowMajorStates, std::size_t samples,
    std::size_t stateDimension, double deltaTime, const Observable& observable,
    const HankelDMDOptions& options = {});

}  // namespace dynsyslab::core
