#include <dynsyslab/core/Integrator.hpp>

#include "IntegratorBackend.hpp"
#include <dynsyslab/core/ControllerCallLog.hpp>
#include <dynsyslab/core/Execution.hpp>

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <limits>
#include <stdexcept>

// One driver for every native ODE method: the output clock (uniform or arbitrary), the
// zero-order-hold control clock, validity and magnitude guards, the per-trajectory step
// budget, checkpoints and pause, and the failure policy. Method references are listed in
// Integrator.hpp and docs/architecture/NATIVE_SOLVER.md.

namespace dynsyslab::core {
namespace detail {

void RhsContext::rhs(double time, const double* state, double* derivative) {
    if (execution) execution->checkpoint();
    if (control) {
        if (control.evaluation == ControlEvaluationMode::sampled) input = heldInput;
        else {
            control.evaluate(trajectory, time, state, input.data());
            ++statistics.controllerCalls;
            if (options.controllerCallLog)
                options.controllerCallLog->record(trajectory, time, state, input.data(),
                    ControllerCallLog::rightHandSide, ControllerCallLog::pending);
        }
        if (control.additiveMatrix.empty()) {
            system.rhsWithInput(time, state, input[0], derivative);
        } else {
            system.rhs(time, state, derivative);
            const std::size_t n = system.descriptor().stateDimension, m = input.size();
            for (std::size_t row = 0; row < n; ++row)
                for (std::size_t column = 0; column < m; ++column)
                    derivative[row] += control.additiveMatrix[row * m + column] * input[column];
        }
    } else {
        system.rhs(time, state, derivative);
    }
    ++statistics.rightHandSideCalls;
}

void RhsContext::accept(double time, const double* state) const {
    if (system.invalid(time, state))
        throw std::runtime_error("native system left its validity domain");
    const std::size_t n = system.descriptor().stateDimension;
    double magnitude = 0.0;
    for (std::size_t index = 0; index < n; ++index) {
        const double value = state[index];
        if (!std::isfinite(value))
            throw std::runtime_error("state guard: non-finite state component "
                                     + std::to_string(index));
        if (index < options.stateLowerBounds.size() && value < options.stateLowerBounds[index])
            throw std::runtime_error("state guard: component " + std::to_string(index)
                                     + " fell below its lower bound");
        if (index < options.stateUpperBounds.size() && value > options.stateUpperBounds[index])
            throw std::runtime_error("state guard: component " + std::to_string(index)
                                     + " exceeded its upper bound");
        magnitude = std::max(magnitude, std::abs(value));
    }
    if (options.maximumStateMagnitude > 0.0 && magnitude > options.maximumStateMagnitude)
        throw std::runtime_error("state guard: state magnitude exceeded its limit");
}

void RhsContext::countStep(bool accepted) {
    if (trajectorySteps >= options.maximumInternalSteps)
        throw std::runtime_error("native integrator exceeded its step budget (max_internal_steps)");
    ++trajectorySteps;
    if (accepted) ++statistics.acceptedSteps; else ++statistics.rejectedSteps;
    // Right-hand-side calls since the previous decision belong to the step decided now.
    if (options.controllerCallLog)
        options.controllerCallLog->resolvePending(accepted ? ControllerCallLog::accepted
                                                           : ControllerCallLog::rejected);
}

}  // namespace detail

namespace {

std::string lower(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char item) { return static_cast<char>(std::tolower(item)); });
    return value;
}

void validate(const ODESystem& system, const std::vector<double>& initialStates,
              std::size_t trajectories, const IntegrationOptions& options) {
    const std::size_t n = system.descriptor().stateDimension;
    if (!trajectories || initialStates.size() != trajectories * n)
        throw std::invalid_argument("initial states do not match the native system dimension");
    if (!std::isfinite(options.initialTime) || !std::isfinite(options.finalTime)
        || options.finalTime <= options.initialTime)
        throw std::invalid_argument("integration requires finite initialTime < finalTime");
    if (options.outputTimes.empty()
        && (!std::isfinite(options.outputStep) || options.outputStep <= 0.0))
        throw std::invalid_argument("integration outputStep must be positive and finite");
    if (!options.outputTimes.empty()) {
        if (options.outputTimes.size() < 2 || options.outputTimes.front() != options.initialTime
            || options.outputTimes.back() != options.finalTime)
            throw std::invalid_argument(
                "integration outputTimes must start at initialTime and end at finalTime");
        for (std::size_t index = 1; index < options.outputTimes.size(); ++index)
            if (!(options.outputTimes[index] > options.outputTimes[index - 1]))
                throw std::invalid_argument("integration outputTimes must be strictly increasing");
    }
    if (!(options.relativeTolerance > 0.0) || !(options.absoluteTolerance > 0.0))
        throw std::invalid_argument("integration tolerances must be positive");
    if (!options.maximumInternalSteps)
        throw std::invalid_argument("integration step budget must be positive");
    if ((!options.stateLowerBounds.empty() && options.stateLowerBounds.size() != n)
        || (!options.stateUpperBounds.empty() && options.stateUpperBounds.size() != n))
        throw std::invalid_argument("state guard bounds must cover every state");
    if (options.checkpointInterval < 0.0 || !std::isfinite(options.checkpointInterval))
        throw std::invalid_argument("checkpoint interval must be nonnegative");
}

void validateControl(const ControlPolicy& control, std::size_t dimension) {
    if (!control) return;
    if (!control.inputDimension) throw std::invalid_argument("control input dimension is zero");
    if (control.vector && control.value)
        throw std::invalid_argument("a control policy has either a scalar or a vector law");
    if (control.additiveMatrix.empty() && control.inputDimension != 1)
        throw std::invalid_argument(
            "a vector input needs additive-state actuation (an n x m additiveMatrix)");
    if (!control.additiveMatrix.empty()
        && control.additiveMatrix.size() != dimension * control.inputDimension)
        throw std::invalid_argument("additive actuation matrix must be state x input");
    if (control.evaluation == ControlEvaluationMode::sampled) {
        if (!std::isfinite(control.sampleTime) || control.sampleTime <= 0.0)
            throw std::invalid_argument("sampled control requires a positive finite sampleTime");
    } else if (control.sampleTime != 0.0) {
        throw std::invalid_argument("continuous-RHS control cannot declare sampleTime");
    }
}

std::size_t checkedProduct(std::size_t left, std::size_t right) {
    if (right && left > (std::numeric_limits<std::size_t>::max)() / right)
        throw ExecutionLimitExceeded("native solution size overflows addressable memory");
    return left * right;
}

// Uniform grid t0 + k h, with the final time appended unless it lies on the grid.
std::vector<double> outputClock(const IntegrationOptions& options) {
    if (!options.outputTimes.empty()) return options.outputTimes;
    const long double duration = static_cast<long double>(options.finalTime) - options.initialTime;
    const long double ratio = duration / options.outputStep;
    if (ratio > 1.0e15L) throw ExecutionLimitExceeded("native solution sample count is not addressable");
    const auto intervals = static_cast<std::size_t>(std::floor(
        ratio + 32.0L * std::numeric_limits<double>::epsilon() * std::max<long double>(1.0L, ratio)));
    std::vector<double> result(intervals + 1);
    for (std::size_t index = 0; index <= intervals; ++index)
        result[index] = std::min(options.initialTime + index * options.outputStep, options.finalTime);
    if (result.back() < options.finalTime - 16.0 * std::numeric_limits<double>::epsilon()
        * std::max(1.0, std::abs(options.finalTime))) result.push_back(options.finalTime);
    else result.back() = options.finalTime;
    return result;
}

double tickTime(double origin, double period, std::size_t tick) {
    const long double value = static_cast<long double>(origin)
        + static_cast<long double>(period) * static_cast<long double>(tick);
    if (!std::isfinite(static_cast<double>(value)))
        throw std::runtime_error("sampled controller clock overflowed");
    return static_cast<double>(value);
}

bool reached(double time, double target) {
    const double tolerance = 16.0 * std::numeric_limits<double>::epsilon()
        * std::max({1.0, std::abs(time), std::abs(target)});
    return time + tolerance >= target;
}

std::unique_ptr<detail::Stepper> makeStepper(const IntegrationOptions& options, std::size_t n) {
    switch (options.method) {
    case IntegrationMethod::radau5: return detail::makeRadauStepper(options, n);
    case IntegrationMethod::bdf: case IntegrationMethod::adams: case IntegrationMethod::lsoda:
        return detail::makeCvodeStepper(options, n);
    default: return detail::makeExplicitStepper(options, n);
    }
}

}  // namespace

IntegrationMethod integrationMethodFromName(const std::string& raw) {
    const auto name = lower(raw);
    if (name == "euler") return IntegrationMethod::euler;
    if (name == "rk4") return IntegrationMethod::rk4;
    if (name == "rk23" || name == "bs23" || name == "ode23") return IntegrationMethod::bogackiShampine23;
    if (name == "rk45" || name == "dopri5" || name == "ode45") return IntegrationMethod::dormandPrince45;
    if (name == "dop853") return IntegrationMethod::dormandPrince853;
    if (name == "radau" || name == "radau5" || name == "radauiia") return IntegrationMethod::radau5;
    if (name == "bdf" || name == "ode15s" || name == "cvode_bdf") return IntegrationMethod::bdf;
    if (name == "adams" || name == "cvode_adams" || name == "ode113") return IntegrationMethod::adams;
    if (name == "lsoda") return IntegrationMethod::lsoda;
    throw std::invalid_argument("unknown ODE method '" + raw
        + "'; native methods: Euler, RK4, RK23, RK45, DOP853, Radau, BDF, Adams, LSODA");
}

std::string integrationMethodName(IntegrationMethod method) {
    switch (method) {
    case IntegrationMethod::euler: return "Euler";
    case IntegrationMethod::rk4: return "RK4";
    case IntegrationMethod::bogackiShampine23: return "RK23";
    case IntegrationMethod::dormandPrince45: return "RK45";
    case IntegrationMethod::dormandPrince853: return "DOP853";
    case IntegrationMethod::radau5: return "Radau";
    case IntegrationMethod::bdf: return "BDF";
    case IntegrationMethod::adams: return "Adams";
    case IntegrationMethod::lsoda: return "LSODA";
    }
    return "unknown";
}

IntegrationFailure::IntegrationFailure(const std::string& reason, NativeSolution partial)
    : std::runtime_error(reason),
      partial_(std::make_shared<const NativeSolution>(std::move(partial))) {}

void ControlPolicy::evaluate(std::size_t trajectory, double time, const double* state,
                             double* input) const {
    if (vector) vector(trajectory, time, state, input);
    else input[0] = value(trajectory, time, state);
}

const double* NativeSolution::state(std::size_t trajectory, std::size_t sample) const {
    if (trajectory >= trajectoryCount || sample >= sampleCount)
        throw std::out_of_range("native solution state index");
    return values.data() + (trajectory * sampleCount + sample) * stateDimension;
}

double NativeSolution::control(std::size_t trajectory, std::size_t sample) const {
    if (controlValues.empty()) throw std::logic_error("native solution has no control history");
    if (trajectory >= trajectoryCount || sample >= sampleCount)
        throw std::out_of_range("native solution control index");
    return controlValues[(trajectory * sampleCount + sample) * inputDimension];
}

TimeSeries NativeSolution::trajectorySeries(std::size_t trajectory) const {
    if (trajectory >= trajectoryCount) throw std::out_of_range("trajectory time series index");
    return TimeSeries::view(time, values.data() + trajectory * sampleCount * stateDimension,
                            sampleCount, stateDimension);
}

TimeSeries NativeSolution::controlSeries(std::size_t trajectory) const {
    if (controlValues.empty()) throw std::logic_error("solution has no control time series");
    if (trajectory >= trajectoryCount) throw std::out_of_range("control time series index");
    return TimeSeries::view(time, controlValues.data() + trajectory * sampleCount * inputDimension,
                            sampleCount, inputDimension);
}

NativeSolution integrate(ODESystem& system, const std::vector<double>& initialStates,
                         std::size_t trajectoryCount, const IntegrationOptions& options,
                         ScalarControlLaw controller, ProgressCallback progress,
                         const ExecutionContext* execution) {
    ControlPolicy policy;
    policy.value = std::move(controller);
    return integrateControlled(system, initialStates, trajectoryCount, options,
                               std::move(policy), std::move(progress), execution);
}

NativeSolution integrateControlled(
    ODESystem& system, const std::vector<double>& initialStates,
    std::size_t trajectoryCount, const IntegrationOptions& options,
    ControlPolicy controller, ProgressCallback progress,
    const ExecutionContext* execution) {
    const std::size_t n = system.descriptor().stateDimension;
    validate(system, initialStates, trajectoryCount, options);
    validateControl(controller, n);
    if (execution) execution->checkpoint();
    const std::size_t m = controller ? controller.inputDimension : 1;
    NativeSolution result;
    result.time = outputClock(options);
    result.trajectoryCount = trajectoryCount;
    result.stateDimension = n;
    result.inputDimension = m;
    result.coordinateSystem = system.descriptor().preferredCoordinates;
    result.sampleCount = result.time.size();
    if (execution && execution->limits.maximumMemoryBytes) {
        const auto rows = checkedProduct(trajectoryCount, result.sampleCount);
        const auto total = checkedProduct(rows, n + (controller ? m : 0)) + result.sampleCount;
        execution->requireMemory(checkedProduct(total, sizeof(double)), "ODE solution buffers");
    }
    result.values.resize(trajectoryCount * result.sampleCount * n);
    if (controller) result.controlValues.resize(trajectoryCount * result.sampleCount * m);
    result.completedSamples.assign(trajectoryCount, result.sampleCount);
    result.terminationReasons.assign(trajectoryCount, {});
    auto stepper = makeStepper(options, n);
    const bool sampled = controller && controller.evaluation == ControlEvaluationMode::sampled;
    auto lastProgress = std::chrono::steady_clock::now();
    double lastFraction = -1.0;
    double lastCheckpoint = options.initialTime;
    const auto store = [&](std::size_t trajectory, std::size_t sample,
                           const std::vector<double>& state, const std::vector<double>& held) {
        std::copy(state.begin(), state.end(),
                  result.values.begin() + (trajectory * result.sampleCount + sample) * n);
        if (!controller) return;
        double* stored = result.controlValues.data() + (trajectory * result.sampleCount + sample) * m;
        if (sampled) std::copy(held.begin(), held.end(), stored);
        else {
            controller.evaluate(trajectory, result.time[sample], state.data(), stored);
            ++result.statistics.controllerCalls;
            if (options.controllerCallLog)
                options.controllerCallLog->record(trajectory, result.time[sample], state.data(), stored,
                    ControllerCallLog::outputSample, ControllerCallLog::accepted);
        }
    };
    const auto logTick = [&](std::size_t trajectory, double time, const std::vector<double>& state,
                             const std::vector<double>& held) {
        if (options.controllerCallLog)
            options.controllerCallLog->record(trajectory, time, state.data(), held.data(),
                ControllerCallLog::sampleTick, ControllerCallLog::accepted);
    };
    const bool fixedStep = options.method == IntegrationMethod::euler
        || options.method == IntegrationMethod::rk4;
    for (std::size_t trajectory = 0; trajectory < trajectoryCount; ++trajectory) {
        if (result.paused) {
            result.completedSamples[trajectory] = 0;
            result.terminationReasons[trajectory] = "paused";
            continue;
        }
        if (execution) execution->checkpoint();
        system.resetRuntimeState();
        std::vector<double> state(initialStates.begin() + trajectory * n,
                                  initialStates.begin() + (trajectory + 1) * n);
        std::vector<double> held(m);
        detail::RhsContext context{system, controller, held, trajectory, result.statistics,
                                   execution, options, 0, std::vector<double>(m)};
        std::size_t completed = 0;
        if (system.invalid(options.initialTime, state.data()))
            throw std::invalid_argument("native initial state is outside the validity domain");
        try {
            if (sampled) {
                controller.evaluate(trajectory, options.initialTime, state.data(), held.data());
                ++result.statistics.controllerCalls;
                logTick(trajectory, options.initialTime, state, held);
            }
            store(trajectory, 0, state, held);
            completed = 1;
            stepper->restart(context, options.initialTime, state);
            std::size_t tick = 1;
            double nextTick = sampled ? tickTime(options.initialTime, controller.sampleTime, tick)
                                      : (std::numeric_limits<double>::infinity)();
            for (std::size_t sample = 1; sample < result.sampleCount; ++sample) {
                double current = result.time[sample - 1];
                const double target = result.time[sample];
                while (!reached(current, target)) {
                    const double end = std::min(target, nextTick);
                    stepper->advance(context, current, end, state);
                    current = end;
                    if (sampled && reached(current, nextTick)) {
                        controller.evaluate(trajectory, nextTick, state.data(), held.data());
                        ++result.statistics.controllerCalls;
                        logTick(trajectory, nextTick, state, held);
                        if (tick == (std::numeric_limits<std::size_t>::max)())
                            throw std::runtime_error("sampled controller tick count overflowed");
                        nextTick = tickTime(options.initialTime, controller.sampleTime, ++tick);
                        stepper->restart(context, current, state);
                    }
                }
                store(trajectory, sample, state, held);
                completed = sample + 1;
                const double fraction = (trajectory * result.sampleCount + sample + 1.0)
                    / (trajectoryCount * result.sampleCount);
                const auto now = std::chrono::steady_clock::now();
                const bool emit = fraction >= 1.0 || fraction - lastFraction >= 0.01
                    || now - lastProgress >= std::chrono::milliseconds(250);
                if (progress && emit)
                    progress({fraction, trajectory, trajectoryCount,
                              result.statistics.acceptedSteps, result.statistics.rejectedSteps});
                if (execution && emit)
                    execution->report({"integrate", "Integrating trajectory", fraction,
                        trajectory * result.sampleCount + sample + 1,
                        trajectoryCount * result.sampleCount});
                if (emit) { lastProgress = now; lastFraction = fraction; }
                if (options.checkpoint && sample + 1 < result.sampleCount
                    && target - lastCheckpoint + 1.0e-12 * std::max(1.0, std::abs(target))
                        >= options.checkpointInterval) {
                    lastCheckpoint = target;
                    const double* controls = controller ? result.controlValues.data()
                        + trajectory * result.sampleCount * m : nullptr;
                    if (!options.checkpoint({trajectory, sample, target, state.data(), n,
                            result.time.data(), result.values.data() + trajectory * result.sampleCount * n,
                            controls, m})) {
                        result.paused = true;
                        break;
                    }
                }
            }
            if (result.paused) {
                result.completedSamples[trajectory] = completed;
                result.terminationReasons[trajectory] = "paused";
            }
        } catch (const ExecutionCancelled&) { throw; }
        catch (const std::exception& failure) {
            if (options.controllerCallLog) options.controllerCallLog->resolvePending(ControllerCallLog::unknown);
            if (!options.truncateOnFailure && !options.keepPartialHistoryOnFailure) throw;
            result.completedSamples[trajectory] = completed;
            result.terminationReasons[trajectory] = failure.what();
            if (!options.truncateOnFailure) {
                // Abort, keeping what was computed: this trajectory up to the failure.
                const double nan = std::numeric_limits<double>::quiet_NaN();
                std::fill(result.values.begin() + (trajectory * result.sampleCount + completed) * n,
                          result.values.end(), nan);
                if (controller)
                    std::fill(result.controlValues.begin()
                              + (trajectory * result.sampleCount + completed) * m,
                              result.controlValues.end(), nan);
                for (std::size_t later = trajectory + 1; later < trajectoryCount; ++later) {
                    result.completedSamples[later] = 0;
                    result.terminationReasons[later] = "not started: an earlier trajectory failed";
                }
                const std::string reason = failure.what();
                throw IntegrationFailure(reason, std::move(result));
            }
        }
        // Calls after the last step decision: a fixed-step method accepts every step.
        if (options.controllerCallLog)
            options.controllerCallLog->resolvePending(fixedStep ? ControllerCallLog::accepted
                                                                : ControllerCallLog::unknown);
        const std::size_t keep = result.completedSamples[trajectory];
        if (keep < result.sampleCount) {
            const double nan = std::numeric_limits<double>::quiet_NaN();
            std::fill(result.values.begin() + (trajectory * result.sampleCount + keep) * n,
                      result.values.begin() + (trajectory + 1) * result.sampleCount * n, nan);
            if (controller)
                std::fill(result.controlValues.begin() + (trajectory * result.sampleCount + keep) * m,
                          result.controlValues.begin() + (trajectory + 1) * result.sampleCount * m, nan);
        }
        lastCheckpoint = options.initialTime;
    }
    return result;
}

}  // namespace dynsyslab::core
