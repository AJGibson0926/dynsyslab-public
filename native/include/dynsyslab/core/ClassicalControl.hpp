#pragma once

#include <dynsyslab/core/LQR.hpp>
#include <dynsyslab/core/PIDCompensator.hpp>
#include <dynsyslab/core/Controllability.hpp>

#include <cstddef>
#include <functional>
#include <limits>
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
    [[nodiscard]] const numerics::ContinuousLqrResult& lastSolution() const noexcept;
    [[nodiscard]] bool stateDependent() const noexcept;

private:
    StateMatrixFunction generator_,input_;
    std::vector<double> statePenalty_,inputPenalty_;
    std::size_t stateDimension_{},inputDimension_{};
    StateTargetFunction target_;
    bool stateDependent_{};
    mutable numerics::ContinuousLqrResult last_;
};

// Finite-horizon LQR: the first-step gain (input x state) and the cost-to-go (state x state).

[[nodiscard]] numerics::FiniteHorizonLqrResult solveFiniteHorizonLQR(
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
    [[nodiscard]] const numerics::FiniteHorizonLqrResult& lastSolution() const noexcept;

private:
    StateMatrixFunction generator_,input_;
    std::vector<double> statePenalty_,inputPenalty_,controlLimits_;
    std::size_t stateDimension_{},inputDimension_{},horizon_{};
    double sampleTime_{};
    StateTargetFunction target_;
    mutable numerics::FiniteHorizonLqrResult last_;
};

// State-feedback PID on e = x_target - x (the "pid" controller kind).
class PIDController {
public:
    PIDController(
        std::vector<double> proportional, std::vector<double> integral,
        std::vector<double> derivative, std::size_t stateDimension,
        std::size_t inputDimension, StateTargetFunction target = {},
        PIDOptions options = {});

    [[nodiscard]] std::vector<double> value(double time, const double* state);
    void reset() noexcept;

private:
    PIDCompensator compensator_;
    std::size_t stateDimension_{};
    StateTargetFunction target_;
};

}  // namespace dynsyslab::core
