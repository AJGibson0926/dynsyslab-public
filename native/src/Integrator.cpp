#include <dynsyslab/core/Integrator.hpp>

#include "IntegratorBackend.hpp"
#include <dynsyslab/core/Execution.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

struct StepContext {
    ODESystem& system;
    const ScalarControlPolicy& controller;
    const double* heldInput{};
    std::size_t trajectory{};
    IntegrationStatistics& statistics;
    const ExecutionContext* execution{};

    void rhs(double time, const double* state, double* derivative) {
        if (execution) execution->checkpoint();
        if (controller) {
            double input{};
            if (controller.evaluation == ControlEvaluationMode::sampled)
                input = *heldInput;
            else {
                input = controller.value(trajectory, time, state);
                ++statistics.controllerCalls;
            }
            system.rhsWithInput(time, state, input, derivative);
        } else {
            system.rhs(time, state, derivative);
        }
        ++statistics.rightHandSideCalls;
    }
};

void validate(const ODESystem& system, const std::vector<double>& initialStates,
              std::size_t trajectories, const IntegrationOptions& options) {
    if (!trajectories || initialStates.size() != trajectories * system.descriptor().stateDimension)
        throw std::invalid_argument("initial states do not match the native system dimension");
    if (!std::isfinite(options.initialTime) || !std::isfinite(options.finalTime)
        || options.finalTime <= options.initialTime)
        throw std::invalid_argument("integration requires finite initialTime < finalTime");
    if (!std::isfinite(options.outputStep) || options.outputStep <= 0.0)
        throw std::invalid_argument("integration outputStep must be positive and finite");
    if (!(options.relativeTolerance > 0.0) || !(options.absoluteTolerance > 0.0))
        throw std::invalid_argument("integration tolerances must be positive");
}

void validateControl(const ScalarControlPolicy& control) {
    if (!control) return;
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

std::size_t checkedSum(std::size_t left, std::size_t right) {
    if (left > (std::numeric_limits<std::size_t>::max)() - right)
        throw ExecutionLimitExceeded("native solution size overflows addressable memory");
    return left + right;
}

std::size_t outputSamples(const IntegrationOptions& options) {
    const long double count = std::ceil(
        static_cast<long double>(options.finalTime - options.initialTime)
        / static_cast<long double>(options.outputStep)) + 1.0L;
    if (count > static_cast<long double>((std::numeric_limits<std::size_t>::max)()))
        throw ExecutionLimitExceeded("native solution sample count is not addressable");
    return static_cast<std::size_t>(count);
}

void requireSolutionMemory(const ExecutionContext* execution,
                           std::size_t trajectories, std::size_t samples,
                           std::size_t dimension, bool controlled) {
    if (!execution || !execution->limits.maximumMemoryBytes) return;
    const auto stateValues = checkedProduct(checkedProduct(trajectories, samples), dimension);
    const auto controlValues = controlled ? checkedProduct(trajectories, samples) : 0;
    const auto totalValues = checkedSum(checkedSum(stateValues, samples), controlValues);
    execution->requireMemory(
        checkedProduct(totalValues, sizeof(double)),
        "ODE solution buffers");
}

void evaluateControl(const ScalarControlLaw& controller, std::size_t trajectory,
                     double time, const double* state, double& output,
                     IntegrationStatistics& statistics) {
    if (!controller) return;
    output = controller(trajectory, time, state);
    ++statistics.controllerCalls;
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

void fixedStep(StepContext& context, IntegrationMethod method, double time, double dt,
               std::vector<double>& state) {
    const std::size_t n = state.size();
    std::vector<double> k1(n), k2(n), k3(n), k4(n), work(n);
    context.rhs(time, state.data(), k1.data());
    if (method == IntegrationMethod::euler) {
        for (std::size_t i = 0; i < n; ++i) state[i] += dt * k1[i];
    } else {
        for (std::size_t i = 0; i < n; ++i) work[i] = state[i] + 0.5*dt*k1[i];
        context.rhs(time + 0.5*dt, work.data(), k2.data());
        for (std::size_t i = 0; i < n; ++i) work[i] = state[i] + 0.5*dt*k2[i];
        context.rhs(time + 0.5*dt, work.data(), k3.data());
        for (std::size_t i = 0; i < n; ++i) work[i] = state[i] + dt*k3[i];
        context.rhs(time + dt, work.data(), k4.data());
        for (std::size_t i = 0; i < n; ++i)
            state[i] += dt * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]) / 6.0;
    }
    ++context.statistics.acceptedSteps;
}

double adaptiveStep(StepContext& context, double time, double dt,
                    std::vector<double>& state, double relativeTolerance,
                    double absoluteTolerance, bool& accepted) {
    const std::size_t n = state.size();
    std::vector<double> k1(n), k2(n), k3(n), k4(n), k5(n), k6(n), k7(n), work(n),
        fifth(n), fourth(n);
    context.rhs(time, state.data(), k1.data());
    for (std::size_t i=0;i<n;++i) work[i]=state[i]+dt*(1.0/5.0)*k1[i];
    context.rhs(time+dt*(1.0/5.0),work.data(),k2.data());
    for (std::size_t i=0;i<n;++i) work[i]=state[i]+dt*(3.0/40.0*k1[i]+9.0/40.0*k2[i]);
    context.rhs(time+dt*(3.0/10.0),work.data(),k3.data());
    for (std::size_t i=0;i<n;++i) work[i]=state[i]+dt*(44.0/45.0*k1[i]-56.0/15.0*k2[i]+32.0/9.0*k3[i]);
    context.rhs(time+dt*(4.0/5.0),work.data(),k4.data());
    for (std::size_t i=0;i<n;++i) work[i]=state[i]+dt*(19372.0/6561.0*k1[i]-25360.0/2187.0*k2[i]+64448.0/6561.0*k3[i]-212.0/729.0*k4[i]);
    context.rhs(time+dt*(8.0/9.0),work.data(),k5.data());
    for (std::size_t i=0;i<n;++i) work[i]=state[i]+dt*(9017.0/3168.0*k1[i]-355.0/33.0*k2[i]+46732.0/5247.0*k3[i]+49.0/176.0*k4[i]-5103.0/18656.0*k5[i]);
    context.rhs(time+dt,work.data(),k6.data());
    for (std::size_t i=0;i<n;++i) fifth[i]=state[i]+dt*(35.0/384.0*k1[i]+500.0/1113.0*k3[i]+125.0/192.0*k4[i]-2187.0/6784.0*k5[i]+11.0/84.0*k6[i]);
    context.rhs(time+dt,fifth.data(),k7.data());
    double error = 0.0;
    for (std::size_t i=0;i<n;++i) {
        fourth[i]=state[i]+dt*(5179.0/57600.0*k1[i]+7571.0/16695.0*k3[i]+393.0/640.0*k4[i]-92097.0/339200.0*k5[i]+187.0/2100.0*k6[i]+1.0/40.0*k7[i]);
        const double scale=absoluteTolerance+relativeTolerance*std::max(std::abs(state[i]),std::abs(fifth[i]));
        error=std::max(error,std::abs(fifth[i]-fourth[i])/scale);
    }
    accepted = error <= 1.0;
    const double factor = error == 0.0 ? 5.0
        : std::clamp(0.9*std::pow(error,-0.2),0.2,5.0);
    if (accepted) state = std::move(fifth);
    return dt * factor;
}

void integrateInterval(StepContext& context, const IntegrationOptions& options,
                       double start, double end, double& candidateStep,
                       std::vector<double>& state) {
    double time = start;
    if (options.method != IntegrationMethod::dormandPrince45) {
        fixedStep(context, options.method, time, end - start, state);
        return;
    }
    while (time < end) {
        if (context.statistics.acceptedSteps + context.statistics.rejectedSteps
            >= options.maximumInternalSteps)
            throw std::runtime_error("native adaptive integrator exceeded maximumInternalSteps");
        const double remaining = end - time;
        const double roundoff = 16.0 * std::numeric_limits<double>::epsilon()
            * std::max({1.0, std::abs(time), std::abs(end)});
        if (remaining <= roundoff) break;
        const double maximum = options.maximumStep > 0.0 ? options.maximumStep : options.outputStep;
        const double dt = std::min({candidateStep, maximum, remaining});
        if (dt <= std::numeric_limits<double>::epsilon() * std::max(1.0, std::abs(time)))
            throw std::runtime_error("native adaptive integrator step underflow");
        bool accepted = false;
        const double suggested = adaptiveStep(context, time, dt, state,
            options.relativeTolerance, options.absoluteTolerance, accepted);
        candidateStep = std::min(maximum, suggested);
        if (accepted) {
            time += dt;
            if (end - time <= roundoff) time = end;
            ++context.statistics.acceptedSteps;
            if (context.system.invalid(time, state.data()))
                throw std::runtime_error("native system left its validity domain");
        } else {
            ++context.statistics.rejectedSteps;
        }
    }
}

}  // namespace

const double* NativeSolution::state(std::size_t trajectory, std::size_t sample) const {
    if (trajectory >= trajectoryCount || sample >= sampleCount)
        throw std::out_of_range("native solution state index");
    return values.data() + (trajectory * sampleCount + sample) * stateDimension;
}

double NativeSolution::control(std::size_t trajectory, std::size_t sample) const {
    if (controlValues.empty()) throw std::logic_error("native solution has no control history");
    if (trajectory >= trajectoryCount || sample >= sampleCount)
        throw std::out_of_range("native solution control index");
    return controlValues[trajectory * sampleCount + sample];
}

TimeSeries NativeSolution::trajectorySeries(std::size_t trajectory) const {
    if (trajectory >= trajectoryCount) throw std::out_of_range("trajectory time series index");
    return TimeSeries::view(time,
        values.data() + trajectory * sampleCount * stateDimension,
        sampleCount, stateDimension);
}

TimeSeries NativeSolution::controlSeries(std::size_t trajectory) const {
    if (controlValues.empty()) throw std::logic_error("solution has no control time series");
    if (trajectory >= trajectoryCount) throw std::out_of_range("control time series index");
    return TimeSeries::view(time, controlValues.data() + trajectory * sampleCount,
                            sampleCount, 1);
}

NativeSolution integrate(ODESystem& system, const std::vector<double>& initialStates,
                         std::size_t trajectoryCount, const IntegrationOptions& options,
                         ScalarControlLaw controller, ProgressCallback progress,
                         const ExecutionContext* execution) {
    return integrateControlled(system, initialStates, trajectoryCount, options,
        {std::move(controller), ControlEvaluationMode::continuousRhs, 0.0},
        std::move(progress), execution);
}

NativeSolution integrateControlled(
    ODESystem& system, const std::vector<double>& initialStates,
    std::size_t trajectoryCount, const IntegrationOptions& options,
    ScalarControlPolicy controller, ProgressCallback progress,
    const ExecutionContext* execution) {
    validate(system, initialStates, trajectoryCount, options);
    validateControl(controller);
    if (execution) execution->checkpoint();
    if (options.method == IntegrationMethod::bdf)
        return detail::integrateBDF(system, initialStates, trajectoryCount, options,
                                    std::move(controller), std::move(progress), execution);
    NativeSolution result;
    result.trajectoryCount = trajectoryCount;
    result.stateDimension = system.descriptor().stateDimension;
    result.coordinateSystem = system.descriptor().preferredCoordinates;
    result.sampleCount = outputSamples(options);
    requireSolutionMemory(execution, trajectoryCount, result.sampleCount,
                          result.stateDimension, static_cast<bool>(controller));
    result.time.resize(result.sampleCount);
    result.values.resize(trajectoryCount * result.sampleCount * result.stateDimension);
    if (controller) result.controlValues.resize(trajectoryCount * result.sampleCount);
    result.completedSamples.assign(trajectoryCount, result.sampleCount);
    result.terminationReasons.assign(trajectoryCount, {});
    for (std::size_t sample = 0; sample < result.sampleCount; ++sample)
        result.time[sample] = std::min(options.initialTime + sample*options.outputStep, options.finalTime);
    auto lastProgress = std::chrono::steady_clock::now();
    double lastFraction = -1.0;
    for (std::size_t trajectory = 0; trajectory < trajectoryCount; ++trajectory) {
        if (execution) execution->checkpoint();
        system.resetRuntimeState();
        std::vector<double> state(initialStates.begin() + trajectory*result.stateDimension,
                                  initialStates.begin() + (trajectory + 1)*result.stateDimension);
        std::copy(state.begin(), state.end(), result.values.begin()
            + trajectory*result.sampleCount*result.stateDimension);
        double heldInput{};
        const bool sampled = controller.evaluation == ControlEvaluationMode::sampled;
        if (controller && sampled)
            evaluateControl(controller.value, trajectory, result.time[0], state.data(),
                            heldInput, result.statistics);
        if (controller) {
            auto& stored = result.controlValues[trajectory * result.sampleCount];
            if (sampled) stored = heldInput;
            else evaluateControl(controller.value, trajectory, result.time[0],
                                 state.data(), stored, result.statistics);
        }
        double candidateStep = options.initialStep > 0.0 ? options.initialStep
            : std::min(options.outputStep, (options.finalTime - options.initialTime) / 100.0);
        StepContext context{system, controller, &heldInput, trajectory,
                            result.statistics, execution};
        std::size_t tick = 1;
        double nextTick = sampled
            ? tickTime(options.initialTime, controller.sampleTime, tick)
            : (std::numeric_limits<double>::infinity)();
        std::size_t completed = 1;
        try {
        for (std::size_t sample = 1; sample < result.sampleCount; ++sample) {
            double current = result.time[sample - 1];
            const double target = result.time[sample];
            while (!reached(current, target)) {
                const double end = std::min(target, nextTick);
                integrateInterval(context, options, current, end, candidateStep, state);
                current = end;
                if (sampled && reached(current, nextTick)) {
                    evaluateControl(controller.value, trajectory, nextTick, state.data(),
                                    heldInput, result.statistics);
                    if (tick == (std::numeric_limits<std::size_t>::max)())
                        throw std::runtime_error("sampled controller tick count overflowed");
                    nextTick = tickTime(options.initialTime, controller.sampleTime, ++tick);
                }
            }
            std::copy(state.begin(), state.end(), result.values.begin()
                + (trajectory*result.sampleCount + sample)*result.stateDimension);
            completed = sample + 1;
            if (controller) {
                auto& stored = result.controlValues[trajectory * result.sampleCount + sample];
                if (sampled) stored = heldInput;
                else evaluateControl(controller.value, trajectory, result.time[sample],
                                     state.data(), stored, result.statistics);
            }
            const double fraction = (trajectory*result.sampleCount + sample + 1.0)
                / (trajectoryCount*result.sampleCount);
            const auto now = std::chrono::steady_clock::now();
            const bool emit = fraction >= 1.0 || fraction - lastFraction >= 0.01
                || now - lastProgress >= std::chrono::milliseconds(250);
            if (progress && emit)
                progress({fraction, trajectory, trajectoryCount, result.statistics.acceptedSteps,
                          result.statistics.rejectedSteps});
            if (execution && emit)
                execution->report({"integrate", "Integrating trajectory", fraction,
                    trajectory * result.sampleCount + sample + 1,
                    trajectoryCount * result.sampleCount});
            if (emit) { lastProgress = now; lastFraction = fraction; }
        }
            } catch (const ExecutionCancelled&) { throw; }
            catch (const std::exception& failure) {
                if (!options.truncateOnFailure) throw;
                const double nan = std::numeric_limits<double>::quiet_NaN();
                for (std::size_t sample = completed; sample < result.sampleCount; ++sample) {
                    std::fill_n(result.values.begin() + (trajectory * result.sampleCount + sample) * result.stateDimension, result.stateDimension, nan);
                    if (controller) result.controlValues[trajectory * result.sampleCount + sample] = nan;
                }
                result.completedSamples[trajectory] = completed;
                result.terminationReasons[trajectory] = failure.what();
            }
    }
    return result;
}

}  // namespace dynsyslab::core
