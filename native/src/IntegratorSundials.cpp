#include "IntegratorBackend.hpp"
#include <dynsyslab/core/Execution.hpp>

#include <cvode/cvode.h>
#include <cvode/cvode_ls.h>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sundials/sundials_context.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <exception>
#include <limits>
#include <stdexcept>
#include <string>

namespace dynsyslab::core::detail {
namespace {

struct CallbackData {
    ODESystem* system{};
    const ScalarControlPolicy* controller{};
    const double* heldInput{};
    std::size_t trajectory{};
    IntegrationStatistics* statistics{};
    const ExecutionContext* execution{};
    std::exception_ptr error;
};

int rightHandSide(sunrealtype time, N_Vector state, N_Vector derivative, void* raw) {
    auto& data = *static_cast<CallbackData*>(raw);
    try {
        if (data.execution) data.execution->checkpoint();
        const double* active = N_VGetArrayPointer(state);
        double* output = N_VGetArrayPointer(derivative);
        if (*data.controller) {
            double input{};
            if (data.controller->evaluation == ControlEvaluationMode::sampled)
                input = *data.heldInput;
            else {
                input = data.controller->value(data.trajectory, time, active);
                ++data.statistics->controllerCalls;
            }
            data.system->rhsWithInput(time, active, input, output);
        } else data.system->rhs(time, active, output);
        ++data.statistics->rightHandSideCalls;
        return 0;
    } catch (...) { data.error = std::current_exception(); return -1; }
}

double tickTime(double origin, double period, std::size_t tick) {
    const long double value = static_cast<long double>(origin)
        + static_cast<long double>(period) * static_cast<long double>(tick);
    if (!std::isfinite(static_cast<double>(value)))
        throw std::runtime_error("sampled BDF controller clock overflowed");
    return static_cast<double>(value);
}

bool reached(double time, double target) {
    const double tolerance = 16.0 * std::numeric_limits<double>::epsilon()
        * std::max({1.0, std::abs(time), std::abs(target)});
    return time + tolerance >= target;
}

void require(int status, const char* operation) {
    if (status < 0) throw std::runtime_error(std::string("SUNDIALS ") + operation
                                             + " failed with code " + std::to_string(status));
}

void accumulateStatistics(void* cvode, IntegrationStatistics& statistics) {
    long int steps = 0, failures = 0;
    require(CVodeGetNumSteps(cvode, &steps), "step statistics");
    require(CVodeGetNumErrTestFails(cvode, &failures), "error statistics");
    statistics.acceptedSteps += static_cast<std::size_t>(std::max(0L, steps));
    statistics.rejectedSteps += static_cast<std::size_t>(std::max(0L, failures));
}

std::vector<double> outputTimes(const IntegrationOptions& options,
                                const ExecutionContext* execution,
                                std::size_t trajectories, std::size_t dimension,
                                bool controlled) {
    const double duration = options.finalTime - options.initialTime;
    const std::size_t intervals = static_cast<std::size_t>(std::floor(
        duration / options.outputStep + 32.0 * std::numeric_limits<double>::epsilon()
        * std::max(1.0, duration / options.outputStep)));
    const std::size_t samples = intervals + 1 + (options.initialTime + intervals
        * options.outputStep < options.finalTime ? 1 : 0);
    if (execution && execution->limits.maximumMemoryBytes) {
        const long double values = static_cast<long double>(samples)
            * (1.0L + static_cast<long double>(trajectories)
               * (dimension + (controlled ? 1.0L : 0.0L)));
        if (values > static_cast<long double>((std::numeric_limits<std::size_t>::max)())
            / sizeof(double))
            throw ExecutionLimitExceeded("native BDF solution size is not addressable");
        execution->requireMemory(static_cast<std::size_t>(values) * sizeof(double),
                                 "BDF solution buffers");
    }
    std::vector<double> result(intervals + 1);
    for (std::size_t index = 0; index <= intervals; ++index)
        result[index] = options.initialTime + index * options.outputStep;
    if (result.back() < options.finalTime - 16.0 * std::numeric_limits<double>::epsilon()
        * std::max(1.0, std::abs(options.finalTime))) result.push_back(options.finalTime);
    else result.back() = options.finalTime;
    return result;
}

struct Resources {
    SUNContext context{}; N_Vector state{}; SUNMatrix matrix{};
    SUNLinearSolver linearSolver{}; void* cvode{};
    ~Resources() {
        if (cvode) CVodeFree(&cvode);
        if (linearSolver) SUNLinSolFree(linearSolver);
        if (matrix) SUNMatDestroy(matrix);
        if (state) N_VDestroy(state);
        if (context) SUNContext_Free(&context);
    }
};

}  // namespace

NativeSolution integrateBDF(ODESystem& system, const std::vector<double>& initialStates,
                            std::size_t trajectories, const IntegrationOptions& options,
                            ScalarControlPolicy controller, ProgressCallback progress,
                            const ExecutionContext* execution) {
    const std::size_t dimension = system.descriptor().stateDimension;
    NativeSolution result; result.time = outputTimes(
        options, execution, trajectories, dimension, static_cast<bool>(controller));
    result.trajectoryCount = trajectories; result.sampleCount = result.time.size();
    result.stateDimension = dimension; result.coordinateSystem = system.descriptor().preferredCoordinates;
    result.values.resize(trajectories * result.sampleCount * dimension);
    if (controller) result.controlValues.resize(trajectories * result.sampleCount);
    result.completedSamples.assign(trajectories, result.sampleCount);
    result.terminationReasons.assign(trajectories, {});
    auto lastProgress = std::chrono::steady_clock::now(); double lastFraction = -1.0;
    for (std::size_t trajectory = 0; trajectory < trajectories; ++trajectory) {
        if (execution) execution->checkpoint();
        system.resetRuntimeState(); Resources resource;
        require(SUNContext_Create(SUN_COMM_NULL, &resource.context), "context creation");
        resource.state = N_VNew_Serial(static_cast<sunindextype>(dimension), resource.context);
        if (!resource.state) throw std::runtime_error("SUNDIALS state allocation failed");
        resource.matrix = SUNDenseMatrix(static_cast<sunindextype>(dimension),
                                         static_cast<sunindextype>(dimension), resource.context);
        if (!resource.matrix) throw std::runtime_error("SUNDIALS matrix allocation failed");
        resource.linearSolver = SUNLinSol_Dense(resource.state, resource.matrix, resource.context);
        if (!resource.linearSolver)
            throw std::runtime_error("SUNDIALS linear-solver allocation failed");
        resource.cvode = CVodeCreate(CV_BDF, resource.context);
        if (!resource.cvode) throw std::runtime_error("SUNDIALS CVODE allocation failed");
        double* active = N_VGetArrayPointer(resource.state);
        std::copy(initialStates.begin() + trajectory * dimension,
                  initialStates.begin() + (trajectory + 1) * dimension, active);
        if (system.invalid(options.initialTime, active))
            throw std::invalid_argument("native BDF initial state is outside the validity domain");
        const bool sampled = controller.evaluation == ControlEvaluationMode::sampled;
        double heldInput{};
        if (controller && sampled) {
            heldInput = controller.value(trajectory, options.initialTime, active);
            ++result.statistics.controllerCalls;
        }
        CallbackData callback{&system, &controller, &heldInput, trajectory,
                              &result.statistics, execution, {}};
        require(CVodeInit(resource.cvode, rightHandSide, options.initialTime, resource.state), "initialization");
        require(CVodeSetUserData(resource.cvode, &callback), "user-data binding");
        require(CVodeSStolerances(resource.cvode, options.relativeTolerance,
                                  options.absoluteTolerance), "tolerance configuration");
        require(CVodeSetLinearSolver(resource.cvode, resource.linearSolver, resource.matrix),
                "dense linear solver binding");
        require(CVodeSetMaxNumSteps(resource.cvode, static_cast<long int>(options.maximumInternalSteps)),
                "maximum-step configuration");
        if (options.initialStep > 0.0)
            require(CVodeSetInitStep(resource.cvode, options.initialStep), "initial-step configuration");
        if (options.maximumStep > 0.0)
            require(CVodeSetMaxStep(resource.cvode, options.maximumStep), "maximum-step configuration");
        std::size_t tick = 1;
        double nextTick = sampled
            ? tickTime(options.initialTime, controller.sampleTime, tick)
            : (std::numeric_limits<double>::infinity)();
        double current = options.initialTime;
        std::size_t completed = 0;
        try {
        for (std::size_t sample = 0; sample < result.sampleCount; ++sample) {
            if (sample) {
                const double target = result.time[sample];
                while (!reached(current, target)) {
                    const double end = std::min(target, nextTick);
                    sunrealtype actual{};
                    const int status = CVode(resource.cvode, end, resource.state,
                                             &actual, CV_NORMAL);
                    if (callback.error) std::rethrow_exception(callback.error);
                    require(status, "integration");
                    current = end;
                    if (system.invalid(actual, active))
                        throw std::runtime_error("native BDF system left its validity domain");
                    if (sampled && reached(current, nextTick)) {
                        accumulateStatistics(resource.cvode, result.statistics);
                        heldInput = controller.value(trajectory, nextTick, active);
                        ++result.statistics.controllerCalls;
                        if (tick == (std::numeric_limits<std::size_t>::max)())
                            throw std::runtime_error("sampled BDF controller tick count overflowed");
                        nextTick = tickTime(options.initialTime, controller.sampleTime, ++tick);
                        require(CVodeReInit(resource.cvode, current, resource.state),
                                "sampled-control reinitialization");
                    }
                }
            }
            std::copy(active, active + dimension, result.values.begin()
                + (trajectory * result.sampleCount + sample) * dimension);
            if (controller) {
                auto& stored = result.controlValues[trajectory * result.sampleCount + sample];
                if (sampled) stored = heldInput;
                else {
                    stored = controller.value(trajectory, result.time[sample], active);
                    ++result.statistics.controllerCalls;
                }
            }
            const double fraction = (trajectory * result.sampleCount + sample + 1.0)
                / (trajectories * result.sampleCount);
            const auto now = std::chrono::steady_clock::now();
            const bool emit = fraction >= 1.0 || fraction-lastFraction >= 0.01
                || now-lastProgress >= std::chrono::milliseconds(250);
            if (progress && emit)
                progress({fraction, trajectory, trajectories, result.statistics.acceptedSteps,
                          result.statistics.rejectedSteps});
            if (execution && emit)
                execution->report({"integrate", "Integrating BDF trajectory", fraction,
                    trajectory * result.sampleCount + sample + 1,
                    trajectories * result.sampleCount});
            if (emit) { lastProgress = now; lastFraction = fraction; }
            completed = sample + 1;
        }
            } catch (const ExecutionCancelled&) { throw; }
            catch (const std::exception& failure) {
                if (!options.truncateOnFailure) throw;
                const double nan = std::numeric_limits<double>::quiet_NaN();
                for (std::size_t sample = completed; sample < result.sampleCount; ++sample) {
                    std::fill_n(result.values.begin() + (trajectory * result.sampleCount + sample) * dimension, dimension, nan);
                    if (controller) result.controlValues[trajectory * result.sampleCount + sample] = nan;
                }
                result.completedSamples[trajectory] = completed;
                result.terminationReasons[trajectory] = failure.what();
            }
        accumulateStatistics(resource.cvode, result.statistics);
    }
    return result;
}

}  // namespace dynsyslab::core::detail
