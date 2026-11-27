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

const numerics::ContinuousLqrResult& StateSpaceLQRController::lastSolution() const noexcept {
    return last_;
}

bool StateSpaceLQRController::stateDependent() const noexcept { return stateDependent_; }

numerics::FiniteHorizonLqrResult solveFiniteHorizonLQR(
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
    return numerics::solveFiniteHorizonLqr(
        generatorValues,inputValues,statePenaltyValues,inputPenaltyValues,
        states,inputs,horizon,sampleTime);
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

const numerics::FiniteHorizonLqrResult& StateSpaceMPCController::lastSolution() const noexcept {
    return last_;
}

PIDCompensator::PIDCompensator(
    std::vector<double> proportional, std::vector<double> integral,
    std::vector<double> derivative, std::size_t errors, std::size_t outputs, PIDOptions options)
    : proportional_(std::move(proportional)),integralGain_(std::move(integral)),
      derivative_(std::move(derivative)),accumulated_(errors),previousError_(errors),
      filteredRate_(errors),errorDimension_(errors),outputDimension_(outputs),
      options_(std::move(options)) {
    dimensions(proportional_, outputs, errors, "PID Kp");
    dimensions(integralGain_, outputs, errors, "PID Ki");
    dimensions(derivative_, outputs, errors, "PID Kd");
    if (!(options_.derivativeFilterTime >= 0.0) || !std::isfinite(options_.derivativeFilterTime))
        throw std::invalid_argument("PID derivative filter time must be nonnegative");
    if (!options_.integratorLimit.empty() && options_.integratorLimit.size() != errors)
        throw std::invalid_argument("PID integrator limit needs one value per error channel");
    for (double limit : options_.integratorLimit)
        if (!(limit > 0.0)) throw std::invalid_argument("PID integrator limits must be positive");
    if (!(options_.outputMinimum < options_.outputMaximum))
        throw std::invalid_argument("PID output limits must satisfy minimum < maximum");
}

std::vector<double> PIDCompensator::combine(const std::vector<double>& error,
    const std::vector<double>& integral, const std::vector<double>& rate) const {
    std::vector<double> output(outputDimension_, 0.0);
    for (std::size_t row = 0; row < outputDimension_; ++row)
        for (std::size_t column = 0; column < errorDimension_; ++column) {
            const std::size_t index = row * errorDimension_ + column;
            output[row] += proportional_[index] * error[column]
                + integralGain_[index] * integral[column] + derivative_[index] * rate[column];
        }
    return output;
}

std::vector<double> PIDCompensator::update(double time, const std::vector<double>& error) {
    if (!std::isfinite(time)) throw std::invalid_argument("PID sample time is invalid");
    if (error.size() != errorDimension_ || std::any_of(error.begin(), error.end(),
        [](double value) { return !std::isfinite(value); }))
        throw std::invalid_argument("PID error has the wrong dimension or is not finite");
    std::vector<double> integral = accumulated_, rate(errorDimension_, 0.0);
    if (initialized_) {
        if (time < previousTime_) throw std::invalid_argument("PID samples must be monotonic");
        const double step = time - previousTime_;
        if (step > 0.0) {
            for (std::size_t index = 0; index < errorDimension_; ++index) {
                integral[index] += step * error[index];
                if (!options_.integratorLimit.empty())
                    integral[index] = std::clamp(integral[index],
                        -options_.integratorLimit[index], options_.integratorLimit[index]);
                const double raw = (error[index] - previousError_[index]) / step;
                // First-order filter, exact for a held raw derivative over the step.
                const double tau = options_.derivativeFilterTime;
                const double blend = tau > 0.0 ? 1.0 - std::exp(-step / tau) : 1.0;
                filteredRate_[index] += blend * (raw - filteredRate_[index]);
            }
        }
        rate = filteredRate_;
    }
    auto output = combine(error, integral, rate);
    if (options_.conditionalIntegration && initialized_) {
        const auto held = combine(error, accumulated_, rate);
        for (std::size_t row = 0; row < outputDimension_; ++row) {
            const bool high = output[row] > options_.outputMaximum && output[row] > held[row];
            const bool low = output[row] < options_.outputMinimum && output[row] < held[row];
            if (high || low) { integral = accumulated_; output = held; break; }
        }
    }
    for (double& value : output)
        value = std::clamp(value, options_.outputMinimum, options_.outputMaximum);
    accumulated_ = std::move(integral); previousError_ = error;
    previousTime_ = time; initialized_ = true;
    return output;
}

void PIDCompensator::reset() noexcept {
    std::fill(accumulated_.begin(), accumulated_.end(), 0.0);
    std::fill(previousError_.begin(), previousError_.end(), 0.0);
    std::fill(filteredRate_.begin(), filteredRate_.end(), 0.0);
    previousTime_ = 0.0; initialized_ = false;
}

PIDController::PIDController(
    std::vector<double> proportional, std::vector<double> integral,
    std::vector<double> derivative, std::size_t states, std::size_t inputs,
    StateTargetFunction target, PIDOptions options)
    : compensator_(std::move(proportional), std::move(integral), std::move(derivative),
                   states, inputs, std::move(options)),
      stateDimension_(states), target_(std::move(target)) {}

std::vector<double> PIDController::value(double time, const double* state) {
    if (!std::isfinite(time)) throw std::invalid_argument("PID sample time is invalid");
    finiteState(state, stateDimension_, "PID");
    const auto target = reference(target_, time, stateDimension_);
    std::vector<double> error(stateDimension_);
    for (std::size_t index = 0; index < stateDimension_; ++index)
        error[index] = target[index] - state[index];
    return compensator_.update(time, error);
}

void PIDController::reset() noexcept { compensator_.reset(); }

}  // namespace dynsyslab::core
