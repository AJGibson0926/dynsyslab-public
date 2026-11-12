#pragma once

#include <dynsyslab/core/KoopmanModel.hpp>
#include <dynsyslab/core/Observable.hpp>

#include <cstddef>
#include <string>
#include <vector>

namespace dynsyslab::core {

struct ControlledData {
    std::vector<double> time;
    std::vector<double> states; // sample x state
    std::vector<double> inputs; // sample x input, or interval x input
    std::size_t samples{};
    std::size_t stateDimension{};
    std::size_t inputDimension{};
};

struct ControlRegressionDiagnostics {
    std::string method;
    std::size_t samples{};
    std::size_t designColumns{};
    std::size_t rank{};
    double conditionNumber{};
    double rootMeanSquareError{};
    double relativeResidual{};
};

struct ConstantControlFit {
    ConstantControlMap representation;
    ControlRegressionDiagnostics diagnostics;
};

struct BilinearControlFit {
    BilinearControlMap representation;
    ControlRegressionDiagnostics diagnostics;
};

[[nodiscard]] ConstantControlFit fitConstantExactControlMap(
    const KoopmanModel& model, const std::vector<double>& states,
    std::size_t samples);

[[nodiscard]] BilinearControlFit fitBilinearExactControlMap(
    const KoopmanModel& model, const std::vector<double>& states,
    std::size_t samples, double ridge = 0.0, bool conjugateComplete = false);

[[nodiscard]] ConstantControlFit fitConstantControlResidual(
    const KoopmanModel& model, const ControlledData& data, double ridge = 0.0);
[[nodiscard]] ConstantControlFit fitConstantControlResidual(
    const KoopmanModel& model, const std::vector<ControlledData>& data,
    double ridge = 0.0);

// Korda--Mezic-style regression: hold A and an observable reconstruction C
// fixed, regress the interval-input convolution in discrete time, then convert
// the fitted input matrix back to the continuous generator convention.
[[nodiscard]] ConstantControlFit fitConstantControlMultistepOutput(
    const KoopmanModel& model, const ControlledData& data,
    const Observable& observable, const ReconstructionResult& reconstruction,
    double ridge = 0.0);
[[nodiscard]] ConstantControlFit fitConstantControlMultistepOutput(
    const KoopmanModel& model, const std::vector<ControlledData>& data,
    const Observable& observable, const ReconstructionResult& reconstruction,
    double ridge = 0.0);

[[nodiscard]] BilinearControlFit fitBilinearControlResidual(
    const KoopmanModel& model, const ControlledData& data,
    double ridge = 0.0, bool conjugateComplete = false);
[[nodiscard]] BilinearControlFit fitBilinearControlResidual(
    const KoopmanModel& model, const std::vector<ControlledData>& data,
    double ridge = 0.0, bool conjugateComplete = false);

struct ControlMapErrorGrid {
    std::vector<double> axis0;
    std::vector<double> axis1;
    std::vector<double> absoluteError; // axis1 x axis0
    std::vector<double> relativeError;
    std::size_t columns{};
    std::size_t rows{};
    double absoluteMean{};
    double absoluteMaximum{};
    double relativeMean{};
    double relativeMaximum{};
};

[[nodiscard]] ControlMapErrorGrid evaluateControlMapError2D(
    const KoopmanModel& model, const std::string& approximation,
    double axis0Minimum, double axis0Maximum,
    double axis1Minimum, double axis1Maximum,
    std::size_t columns = 101, std::size_t rows = 101,
    const std::vector<double>& fixedState = {});

}  // namespace dynsyslab::core
