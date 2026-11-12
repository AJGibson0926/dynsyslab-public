#include <dynsyslab/core/ClassicalControl.hpp>
#include <dynsyslab/numerics/Control.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

namespace dynsyslab::core {
namespace {

void dimensions(const std::vector<double>& matrix, std::size_t rows,
                std::size_t columns, const char* name) {
    if (!rows || !columns || matrix.size() != rows * columns)
        throw std::invalid_argument(std::string(name) + " dimensions are inconsistent");
    if (std::any_of(matrix.begin(),matrix.end(),
        [](double value){return !std::isfinite(value);}))
        throw std::invalid_argument(std::string(name) + " must be finite");
}

std::vector<double> reference(
    const StateTargetFunction& target, double time, std::size_t dimension) {
    std::vector<double> result(dimension, 0.0);
    if (target) result = target(time);
    if (result.size() != dimension || std::any_of(result.begin(), result.end(),
        [](double value) { return !std::isfinite(value); }))
        throw std::invalid_argument("controller target dimension or values are invalid");
    return result;
}

std::vector<double> feedback(const std::vector<double>& gain,
                             std::size_t inputs, std::size_t states,
                             const double* state,
                             const std::vector<double>& target) {
    std::vector<double> result(inputs, 0.0);
    for (std::size_t input = 0; input < inputs; ++input)
        for (std::size_t index = 0; index < states; ++index)
            result[input] -= gain[input * states + index]
                * (state[index] - target[index]);
    return result;
}

void limits(std::vector<double>& value, const std::vector<double>& requested) {
    if (requested.empty()) return;
    for (std::size_t index = 0; index < value.size(); ++index) {
        const double limit = requested.size() == 1 ? requested.front() : requested[index];
        value[index] = std::clamp(value[index], -limit, limit);
    }
}

void finiteState(const double* state, std::size_t dimension, const char* name) {
    if (!state) throw std::invalid_argument(std::string(name) + " state is null");
    for (std::size_t index = 0; index < dimension; ++index)
        if (!std::isfinite(state[index]))
            throw std::invalid_argument(std::string(name) + " state must be finite");
}

}  // namespace

StateSpaceLQRController::StateSpaceLQRController(
    std::vector<double> generator, std::vector<double> input,
    std::vector<double> statePenalty, std::vector<double> inputPenalty,
    std::size_t states, std::size_t inputs, StateTargetFunction target)
    : StateSpaceLQRController(
        [value = std::move(generator)](const double*) { return value; },
        [value = std::move(input)](const double*) { return value; },
        std::move(statePenalty), std::move(inputPenalty), states, inputs,
        std::move(target)) {
    stateDependent_ = false;
}

StateSpaceLQRController::StateSpaceLQRController(
    StateMatrixFunction generator, StateMatrixFunction input,
    std::vector<double> statePenalty, std::vector<double> inputPenalty,
    std::size_t states, std::size_t inputs, StateTargetFunction target)
    : generator_(std::move(generator)),input_(std::move(input)),
      statePenalty_(std::move(statePenalty)),inputPenalty_(std::move(inputPenalty)),
      stateDimension_(states),inputDimension_(inputs),target_(std::move(target)),
      stateDependent_(true) {
    if (!generator_ || !input_) throw std::invalid_argument("LQR matrix callback is empty");
    dimensions(statePenalty_, states, states, "LQR Q");
    dimensions(inputPenalty_, inputs, inputs, "LQR R");
}

std::vector<double> StateSpaceLQRController::value(
    double time, const double* state) const {
    finiteState(state, stateDimension_, "LQR");
    if (stateDependent_ || last_.gain.empty()) {
        auto generator = generator_(state), input = input_(state);
        last_ = solveContinuousLQR(generator, input, statePenalty_, inputPenalty_,
                                   stateDimension_, inputDimension_);
    }
    return feedback(last_.gain, inputDimension_, stateDimension_, state,
                    reference(target_, time, stateDimension_));
}

const ContinuousLQRSolution& StateSpaceLQRController::lastSolution() const noexcept {
    return last_;
}

bool StateSpaceLQRController::stateDependent() const noexcept { return stateDependent_; }

FiniteHorizonLQRSolution solveFiniteHorizonLQR(
    const std::vector<double>& generatorValues, const std::vector<double>& inputValues,
    const std::vector<double>& statePenaltyValues,
    const std::vector<double>& inputPenaltyValues, std::size_t states,
    std::size_t inputs, std::size_t horizon, double sampleTime) {
    dimensions(generatorValues, states, states, "MPC A");
    dimensions(inputValues, states, inputs, "MPC B");
    dimensions(statePenaltyValues, states, states, "MPC Q");
    dimensions(inputPenaltyValues, inputs, inputs, "MPC R");
    if (!horizon || !std::isfinite(sampleTime) || sampleTime <= 0.0)
        throw std::invalid_argument("MPC horizon and sample time must be positive");
    auto solved=numerics::solveFiniteHorizonLqr(
        generatorValues,inputValues,statePenaltyValues,inputPenaltyValues,
        states,inputs,horizon,sampleTime);
    FiniteHorizonLQRSolution result;
    result.stateDimension = states; result.inputDimension = inputs;
    result.gain=std::move(solved.gain);
    result.cost=std::move(solved.cost);
    return result;
}

StateSpaceMPCController::StateSpaceMPCController(
    StateMatrixFunction generator, StateMatrixFunction input,
    std::vector<double> statePenalty, std::vector<double> inputPenalty,
    std::size_t states, std::size_t inputs, StateTargetFunction target,
    std::size_t horizon, double sampleTime, std::vector<double> controlLimits)
    : generator_(std::move(generator)),input_(std::move(input)),
      statePenalty_(std::move(statePenalty)),inputPenalty_(std::move(inputPenalty)),
      controlLimits_(std::move(controlLimits)),stateDimension_(states),
      inputDimension_(inputs),horizon_(horizon),sampleTime_(sampleTime),
      target_(std::move(target)) {
    if (!generator_ || !input_) throw std::invalid_argument("MPC matrix callback is empty");
    dimensions(statePenalty_, states, states, "MPC Q");
    dimensions(inputPenalty_, inputs, inputs, "MPC R");
    if (!controlLimits_.empty() && controlLimits_.size() != 1
        && controlLimits_.size() != inputs)
        throw std::invalid_argument("MPC limits require one or one-per-input values");
    for (const auto limit : controlLimits_)
        if (!std::isfinite(limit) || limit <= 0.0)
            throw std::invalid_argument("MPC limits must be positive and finite");
    if (!horizon_ || !std::isfinite(sampleTime_) || sampleTime_ <= 0.0)
        throw std::invalid_argument("MPC horizon and sample time must be positive");
}

std::vector<double> StateSpaceMPCController::value(
    double time, const double* state) const {
    finiteState(state, stateDimension_, "MPC");
    last_ = solveFiniteHorizonLQR(generator_(state), input_(state), statePenalty_,
        inputPenalty_, stateDimension_, inputDimension_, horizon_, sampleTime_);
    auto result = feedback(last_.gain, inputDimension_, stateDimension_, state,
                           reference(target_, time, stateDimension_));
    limits(result, controlLimits_);
    return result;
}

const FiniteHorizonLQRSolution& StateSpaceMPCController::lastSolution() const noexcept {
    return last_;
}

PIDController::PIDController(
    std::vector<double> proportional, std::vector<double> integral,
    std::vector<double> derivative, std::size_t states, std::size_t inputs,
    StateTargetFunction target)
    : proportional_(std::move(proportional)),integralGain_(std::move(integral)),
      derivative_(std::move(derivative)),accumulated_(states),previousError_(states),
      stateDimension_(states),inputDimension_(inputs),target_(std::move(target)) {
    dimensions(proportional_, inputs, states, "PID Kp");
    dimensions(integralGain_, inputs, states, "PID Ki");
    dimensions(derivative_, inputs, states, "PID Kd");
}

std::vector<double> PIDController::value(double time, const double* state) {
    if (!std::isfinite(time)) throw std::invalid_argument("PID sample time is invalid");
    finiteState(state, stateDimension_, "PID");
    const auto target = reference(target_, time, stateDimension_);
    std::vector<double> error(stateDimension_), rate(stateDimension_);
    for (std::size_t index = 0; index < stateDimension_; ++index)
        error[index] = target[index] - state[index];
    if (initialized_) {
        if (time < previousTime_) throw std::invalid_argument("PID samples must be monotonic");
        const double step = time - previousTime_;
        if (step > 0.0)
            for (std::size_t index = 0; index < stateDimension_; ++index) {
                accumulated_[index] += step * error[index];
                rate[index] = (error[index] - previousError_[index]) / step;
            }
    }
    previousError_ = error; previousTime_ = time; initialized_ = true;
    std::vector<double> output(inputDimension_, 0.0);
    for (std::size_t input = 0; input < inputDimension_; ++input)
        for (std::size_t stateIndex = 0; stateIndex < stateDimension_; ++stateIndex) {
            const std::size_t index = input * stateDimension_ + stateIndex;
            output[input] += proportional_[index] * error[stateIndex]
                + integralGain_[index] * accumulated_[stateIndex]
                + derivative_[index] * rate[stateIndex];
        }
    return output;
}

void PIDController::reset() noexcept {
    std::fill(accumulated_.begin(), accumulated_.end(), 0.0);
    std::fill(previousError_.begin(), previousError_.end(), 0.0);
    previousTime_ = 0.0; initialized_ = false;
}

}  // namespace dynsyslab::core
