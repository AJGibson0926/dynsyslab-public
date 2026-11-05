#pragma once

#include <dynsyslab/core/LQR.hpp>

#include <cstddef>
#include <functional>
#include <vector>

namespace dynsyslab::core {

using StateMatrixFunction = std::function<std::vector<double>(const double*)>;

class StateSpaceLQRController {
public:
    StateSpaceLQRController(
        std::vector<double> generator, std::vector<double> input,
        std::vector<double> statePenalty, std::vector<double> inputPenalty,
        std::size_t stateDimension, std::size_t inputDimension,
        StateTargetFunction target = {});
    StateSpaceLQRController(
        StateMatrixFunction generator, StateMatrixFunction input,
        std::vector<double> statePenalty, std::vector<double> inputPenalty,
        std::size_t stateDimension, std::size_t inputDimension,
        StateTargetFunction target = {});

    // Returns the input u. The dynamical system applies its own input map B(x).
    [[nodiscard]] std::vector<double> value(double time, const double* state) const;
    [[nodiscard]] const ContinuousLQRSolution& lastSolution() const noexcept;
    [[nodiscard]] bool stateDependent() const noexcept;

private:
    StateMatrixFunction generator_,input_;
    std::vector<double> statePenalty_,inputPenalty_;
    std::size_t stateDimension_{},inputDimension_{};
    StateTargetFunction target_;
    bool stateDependent_{};
    mutable ContinuousLQRSolution last_;
};

struct FiniteHorizonLQRSolution {
    std::vector<double> gain; // input x state, row major
    std::vector<double> cost; // state x state, row major
    std::size_t stateDimension{};
    std::size_t inputDimension{};
};

[[nodiscard]] FiniteHorizonLQRSolution solveFiniteHorizonLQR(
    const std::vector<double>& generator, const std::vector<double>& input,
    const std::vector<double>& statePenalty, const std::vector<double>& inputPenalty,
    std::size_t stateDimension, std::size_t inputDimension,
    std::size_t horizon, double sampleTime);

class StateSpaceMPCController {
public:
    StateSpaceMPCController(
        StateMatrixFunction generator, StateMatrixFunction input,
        std::vector<double> statePenalty, std::vector<double> inputPenalty,
        std::size_t stateDimension, std::size_t inputDimension,
        StateTargetFunction target, std::size_t horizon = 20,
        double sampleTime = 1.0e-3, std::vector<double> controlLimits = {});

    [[nodiscard]] std::vector<double> value(double time, const double* state) const;
    [[nodiscard]] const FiniteHorizonLQRSolution& lastSolution() const noexcept;

private:
    StateMatrixFunction generator_,input_;
    std::vector<double> statePenalty_,inputPenalty_,controlLimits_;
    std::size_t stateDimension_{},inputDimension_{},horizon_{};
    double sampleTime_{};
    StateTargetFunction target_;
    mutable FiniteHorizonLQRSolution last_;
};

class PIDController {
public:
    PIDController(
        std::vector<double> proportional, std::vector<double> integral,
        std::vector<double> derivative, std::size_t stateDimension,
        std::size_t inputDimension, StateTargetFunction target = {});

    [[nodiscard]] std::vector<double> value(double time, const double* state);
    void reset() noexcept;

private:
    std::vector<double> proportional_,integralGain_,derivative_;
    std::vector<double> accumulated_,previousError_;
    std::size_t stateDimension_{},inputDimension_{};
    StateTargetFunction target_;
    double previousTime_{};
    bool initialized_{};
};

}  // namespace dynsyslab::core
