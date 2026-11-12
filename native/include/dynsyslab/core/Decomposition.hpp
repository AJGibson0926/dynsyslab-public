#pragma once

#include <complex>
#include <cstddef>
#include <string>
#include <vector>

#include <dynsyslab/core/Observable.hpp>

namespace dynsyslab::core {

enum class DictionaryKind { polynomial, radialBasis };

struct DictionaryOptions {
    DictionaryKind kind{DictionaryKind::polynomial};
    std::size_t polynomialDegree{2};
    bool includeConstant{true};
    std::vector<double> radialBasisCenters; // center x state
    std::size_t centerCount{};
    double radialBasisWidth{1.0};
};

struct EDMDOptions {
    DictionaryOptions dictionary;
    double regularization{};
};

struct EDMDResult {
    std::vector<std::complex<double>> discreteEigenvalues;
    std::vector<std::complex<double>> continuousEigenvalues;
    std::vector<std::complex<double>> modes; // state x mode
    std::vector<std::complex<double>> eigenfunctionValues; // sample x mode
    std::vector<std::complex<double>> eigenfunctionCoefficients; // feature x mode
    std::vector<double> operatorMatrix; // feature x feature
    std::size_t samples{};
    std::size_t stateDimension{};
    std::size_t featureCount{};
    double deltaTime{};
};

[[nodiscard]] EDMDResult extendedDynamicModeDecomposition(
    const std::vector<double>& states, std::size_t samples,
    std::size_t stateDimension, double deltaTime,
    const EDMDOptions& options = {});

struct HAVOKOptions {
    std::size_t delays{100};
    std::size_t rank{15};
    std::size_t observableIndex{};
};

struct HAVOKResult {
    std::vector<double> generator; // (rank-1) x (rank-1)
    std::vector<double> input;     // (rank-1) x 1
    std::vector<double> forcing;   // embedded samples
    std::vector<double> coordinates; // embedded samples x rank
    std::vector<double> coordinateTime;
    std::vector<double> leftSingularVectors; // delays x rank
    std::vector<double> singularValues;
    std::vector<double> rightSingularVectors; // rank x embedded samples
    std::vector<double> observable;
    std::size_t delays{};
    std::size_t rank{};
    std::size_t embeddedSamples{};
    double deltaTime{};
};

[[nodiscard]] HAVOKResult havokDecomposition(
    const std::vector<double>& states, const std::vector<double>& time,
    std::size_t samples, std::size_t stateDimension,
    const HAVOKOptions& options = {});

[[nodiscard]] HAVOKResult havokDecomposition(
    const std::vector<double>& states, const std::vector<double>& time,
    std::size_t samples, std::size_t stateDimension,
    const Observable& observable, const HAVOKOptions& options = {});

}  // namespace dynsyslab::core
