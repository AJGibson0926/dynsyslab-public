#include "MeasurementMap.hpp"
#include <dynsyslab/workflow/ControlFactory.hpp>
#include "KlqrAugmentation.hpp"
#include "ModelBindingOptions.hpp"
#include "BalancedTruncationSettings.hpp"
#include "BilinearCertificateBinding.hpp"
#include "HistoryControl.hpp"
#include "MeasuredJetControl.hpp"
#include "IoDelayControl.hpp"
#include "DerivativeSettings.hpp"
#include "EstimatedControl.hpp"
#include <dynsyslab/core/HistoryObservationArchive.hpp>

#include <dynsyslab/core/ClassicalControl.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/Text.hpp>
#include <dynsyslab/workflow/ControlTargets.hpp>

#include "PrescribedControl.hpp"
#include "ControlFactoryInternal.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

namespace dynsyslab::workflow {
namespace control_factory_detail {

std::vector<double> numbers(const json::Value& value) {
    std::vector<double> result;
    for (const auto& item : value.array()) result.push_back(item.number());
    return result;
}

std::vector<double> squareMatrix(const json::Value& source,
                                 std::size_t dimension, const char* name) {
    const auto values = numbers(source);
    std::vector<double> result(dimension * dimension, 0.0);
    if (values.empty()) {
        const double scalar = source.number(1.0);
        for (std::size_t index = 0; index < dimension; ++index)
            result[index * dimension + index] = scalar;
    } else if (values.size() == dimension) {
        for (std::size_t index = 0; index < dimension; ++index)
            result[index * dimension + index] = values[index];
    } else if (values.size() == dimension * dimension) {
        result = values;
    } else {
        throw std::invalid_argument(std::string(name)
            + " must be scalar, diagonal, or square");
    }
    if (std::any_of(result.begin(), result.end(),
        [](double value) { return !std::isfinite(value); }))
        throw std::invalid_argument(std::string(name) + " must be finite");
    return result;
}

std::size_t positiveCount(const json::Value& source, std::size_t fallback,
                          const char* name) {
    const double value = source.number(static_cast<double>(fallback));
    if (!std::isfinite(value) || value < 1.0 || std::floor(value) != value
        || value > 9007199254740991.0)
        throw std::invalid_argument(std::string(name) + " must be a positive integer");
    return static_cast<std::size_t>(value);
}


Clock resolveClock(const std::string& kind, const json::Value& solver,
                   const json::Value& settings, double fallback,
                   bool sampledRequired) {
    const std::string requested = solver.at("control_evaluation").string("automatic");
    const bool sampled = requested == "sampled"
        || (requested == "automatic" && sampledRequired);
    if (requested != "automatic" && requested != "sampled"
        && requested != "continuous_rhs")
        throw std::invalid_argument("control_evaluation must be automatic, sampled, or continuous_rhs");
    if (requested == "continuous_rhs" && sampledRequired)
        throw std::invalid_argument(kind
            + " requires a sampled controller clock; continuous RHS evaluation is unsafe");
    Clock result;
    if (!sampled) return result;
    result.mode = dynsyslab::core::ControlEvaluationMode::sampled;
    result.label = "sampled";
    result.sampleTime = !solver.at("control_sample_time").isNull()
        ? solver.at("control_sample_time").number()
        : (!settings.at("sample_time").isNull()
            ? settings.at("sample_time").number() : fallback);
    if (!std::isfinite(result.sampleTime) || result.sampleTime <= 0.0)
        throw std::invalid_argument("sampled control requires a positive control sample time");
    return result;
}

}  // namespace control_factory_detail
using namespace control_factory_detail;
namespace {

using dynsyslab::core::ScalarControlLaw;
using dynsyslab::core::StateTargetFunction;

std::vector<double> explicitMatrix(const json::Value& source,
                                   std::size_t rows, std::size_t columns,
                                   const char* name) {
    auto result = numbers(source);
    if (result.size() != rows * columns)
        throw std::invalid_argument(std::string(name)
            + " must contain exactly " + std::to_string(rows * columns)
            + " row-major values");
    if (std::any_of(result.begin(), result.end(),
        [](double value) { return !std::isfinite(value); }))
        throw std::invalid_argument(std::string(name) + " must be finite");
    return result;
}

std::vector<double> gainMatrix(const json::Value& source, std::size_t inputs,
                               std::size_t states, double fallback,
                               const char* name) {
    auto values = numbers(source);
    if (values.empty()) {
        values.assign(inputs * states, 0.0);
        const double scalar = source.number(fallback);
        for (std::size_t index = 0; index < std::min(inputs, states); ++index)
            values[index * states + index] = scalar;
    } else if (inputs == 1 && values.size() == states) {
        // A state-length row is the natural scalar-input PID form.
    } else if (states == 1 && values.size() == inputs) {
        // An input-length column is the natural scalar-state form.
    } else if (values.size() != inputs * states) {
        throw std::invalid_argument(std::string(name)
            + " must be scalar or an input-by-state matrix");
    }
    if (std::any_of(values.begin(), values.end(),
        [](double value) { return !std::isfinite(value); }))
        throw std::invalid_argument(std::string(name) + " must be finite");
    return values;
}

std::vector<double> activeWindow(const json::Value& settings) {
    auto active = numbers(settings.at("active_time_window"));
    if (!active.empty() && (active.size() != 2 || active[0] > active[1]
        || !std::isfinite(active[0]) || !std::isfinite(active[1])))
        throw std::invalid_argument("active_time_window must be [start, stop]");
    return active;
}

ScalarControlLaw withWindow(ScalarControlLaw law, std::vector<double> active) {
    if (active.empty()) return law;
    return [law = std::move(law), active = std::move(active)](
        std::size_t trajectory, double time, const double* state) {
        if (time < active[0] || time > active[1]) return 0.0;
        return law(trajectory, time, state);
    };
}

void requireScalarInput(const dynsyslab::core::ODESystem& system,
                        const double* state) {
    const auto matrix = system.inputMatrix(state);
    if (matrix.size() != system.descriptor().stateDimension)
        throw std::invalid_argument(
            "current native ODE integrator supports one scalar system input");
}

std::vector<double> linearizationPoint(const json::Value& settings,
                                       const std::vector<double>& fallback,
                                       std::size_t dimension) {
    auto point = numbers(settings.at("linearization_state"));
    if (point.empty()) point = fallback;
    if (point.size() < dimension)
        throw std::invalid_argument("linearization_state dimension mismatch");
    point.resize(dimension);
    if (std::any_of(point.begin(), point.end(),
        [](double value) { return !std::isfinite(value); }))
        throw std::invalid_argument("linearization_state must be finite");
    return point;
}

template<class Controller, class Builder>
ScalarControlLaw controllerBank(std::size_t count, Builder builder) {
    auto bank = std::make_shared<std::vector<std::shared_ptr<Controller>>>();
    bank->reserve(count);
    for (std::size_t index = 0; index < count; ++index)
        bank->push_back(builder());
    return [bank = std::move(bank)](std::size_t trajectory, double time,
                                    const double* state) {
        if (trajectory >= bank->size())
            throw std::out_of_range("controller trajectory index is invalid");
        return bank->at(trajectory)->value(time, state).front();
    };
}

}  // namespace

ControllerSetup buildController(
    const json::Value& source, const json::Value& solver,
    const dynsyslab::core::ODESystem& system, std::size_t trajectoryCount,
    double defaultSampleTime, const std::vector<double>& fallbackState) {
    const std::string kind = source.at("kind").string("none");
    if (kind == "none") return {};
    if (!trajectoryCount) throw std::invalid_argument("controller needs a trajectory");
    const auto& settings = source.at("settings");
    if(!settings.at("estimator").isNull()) {
        // applied_controller: the KLQE -> KLQR runs in shadow beside a reference drive that the
        // plant actually receives (any non-estimating controller, evaluated on the true state).
        std::function<double(std::size_t,double,const double*)> applied;
        if(const auto& spec=settings.at("applied_controller");!spec.isNull()) {
            if(!spec.at("settings").at("estimator").isNull()||!spec.at("settings").at("applied_controller").isNull())
                throw std::invalid_argument("applied_controller must not estimate or shadow");
            applied=buildController(spec,solver,system,trajectoryCount,defaultSampleTime,fallbackState).policy.value;
            if(!applied)throw std::invalid_argument("applied_controller needs an input policy");
        }
        return detail::buildEstimatedController(source,solver,system,trajectoryCount,defaultSampleTime,
                                                std::move(applied));
    }
    const std::size_t states = system.descriptor().stateDimension;
    const bool prescribedKind = detail::isPrescribedControl(kind);
    const auto& targetSpec=source.at("target");
    const auto targetSpace=targetSpec.at("space").string("state");
    if(targetSpace!="state"&&targetSpace!="observable"&&targetSpace!="lifted")
        throw std::invalid_argument("control target space must be state, observable or lifted");
    const bool observableTarget=targetSpace=="observable";
    const bool liftedTarget=targetSpace=="lifted";  // z_ref in the K-ROM's real split coordinates
    if(observableTarget&&kind!="klqr"&&kind!="kmpc"&&kind!="io_delay_lqr")
        throw std::invalid_argument("observable-coordinate targets require Koopman control");
    if(liftedTarget&&kind!="klqr"&&kind!="kmpc")
        throw std::invalid_argument("lifted targets require KLQR or KMPC");
    StateTargetFunction target = prescribedKind || observableTarget || liftedTarget ? StateTargetFunction{}
        : resolveControlTarget(targetSpec, states);
    const auto point = linearizationPoint(settings, fallbackState, states);
    // actuation "additive_state": the input enters as xdot = f(x) + B u with B the n x 1
    // input_matrix (or the explicit LQR B), instead of through the system's input channel.
    const auto actuation = settings.at("actuation").string("system_input");
    if (actuation != "system_input" && actuation != "additive_state")
        throw std::invalid_argument("controller actuation must be system_input or additive_state");
    const bool additive = actuation == "additive_state";
    if (additive && !prescribedKind && kind != "pid"
        && !((kind == "lqr" || kind == "mpc" || kind == "klqr" || kind == "kmpc")
             && !(settings.at("B").isNull() && settings.at("input_matrix").isNull())))
        throw std::invalid_argument("additive_state actuation supports prescribed inputs, pid, and lqr, mpc, klqr or kmpc with an explicit B");
    if (!additive) requireScalarInput(system, point.data());

    ControllerSetup result;
    bool historyClock=false;
    if (!prescribedKind && buildMultiInputController(kind, settings, solver, system, trajectoryCount,
                                                     defaultSampleTime, target, result))
        return result;
    if (prescribedKind) {
        auto prescribed = detail::buildPrescribedControl(kind, settings, system);
        result.kind = kind; result.family = "prescribed";
        result.trajectoryMetadata = std::move(prescribed.trajectoryMetadata);
        result.policy.value = withWindow(std::move(prescribed.law),
                                         std::move(prescribed.activeWindow));
        result.settingsCoordinates = std::move(prescribed.settingsCoordinates);
    } else if (kind == "lqr" || kind == "sdre") {
        const auto Q = squareMatrix(settings.at("Q"), states, "Q");
        const auto R = squareMatrix(settings.at("R"), 1, "R");
        const auto generator = [&system, states](const double* state) {
            auto value = system.jacobian(0.0, state);
            if (value.size() != states * states)
                throw std::invalid_argument("system Jacobian dimension mismatch");
            return value;
        };
        const auto input = [&system, states](const double* state) {
            auto value = system.inputMatrix(state);
            if (value.size() != states)
                throw std::invalid_argument("system input matrix dimension mismatch");
            return value;
        };
        if (kind == "lqr") {
            const bool hasA = !settings.at("A").isNull();
            const bool hasB = !settings.at("B").isNull();
            if (hasA != hasB)
                throw std::invalid_argument(
                    "classical LQR requires both explicit A and B, or neither");
            const auto A = hasA ? explicitMatrix(settings.at("A"), states, states, "A")
                                : generator(point.data());
            const auto B = hasB ? explicitMatrix(settings.at("B"), states, 1, "B")
                                : input(point.data());
            result.policy.value = controllerBank<dynsyslab::core::StateSpaceLQRController>(
                trajectoryCount, [A, B, Q, R, states, target] {
                    return std::make_shared<dynsyslab::core::StateSpaceLQRController>(
                        A, B, Q, R, states, 1, target);
                });
        } else if (!settings.at("sdc_A").isNull()) {
            // User SDC factorisation f(x) = A(x) x, B(x) (ControlFactoryClassical.cpp).
            const auto userA = sdcFactor(settings.at("sdc_A"), states, states, system, "sdc_A");
            const auto userB = settings.at("sdc_B").isNull() ? dynsyslab::core::StateMatrixFunction(input)
                : sdcFactor(settings.at("sdc_B"), states, 1, system, "sdc_B");
            result.policy.value = controllerBank<dynsyslab::core::StateSpaceLQRController>(
                trajectoryCount, [userA, userB, Q, R, states, target] {
                    return std::make_shared<dynsyslab::core::StateSpaceLQRController>(
                        userA, userB, Q, R, states, 1, target);
                });
        } else {
            result.policy.value = controllerBank<dynsyslab::core::StateSpaceLQRController>(
                trajectoryCount, [generator, input, Q, R, states, target] {
                    return std::make_shared<dynsyslab::core::StateSpaceLQRController>(
                        generator, input, Q, R, states, 1, target);
                });
        }
        result.kind = kind; result.family = "state_space";
    } else if (kind == "mpc") {
        const auto Q = squareMatrix(settings.at("Q"), states, "Q");
        const auto R = squareMatrix(settings.at("R"), 1, "R");
        auto limits = numbers(settings.at("control_limit"));
        if (limits.empty() && !settings.at("control_limit").isNull())
            limits.push_back(settings.at("control_limit").number());
        const auto generator = [&system, states](const double* state) {
            auto value = system.jacobian(0.0, state);
            if (value.size() != states * states)
                throw std::invalid_argument("system Jacobian dimension mismatch");
            return value;
        };
        const auto input = [&system, states](const double* state) {
            auto value = system.inputMatrix(state);
            if (value.size() != states)
                throw std::invalid_argument("system input matrix dimension mismatch");
            return value;
        };
        const auto clock = resolveClock(kind, solver, settings, defaultSampleTime, true);
        result.policy.value = controllerBank<dynsyslab::core::StateSpaceMPCController>(
            trajectoryCount, [generator, input, Q, R, states, target, limits,
                              horizon = positiveCount(
                                  settings.at("horizon"), 20, "MPC horizon"),
                              sampleTime = clock.sampleTime] {
                return std::make_shared<dynsyslab::core::StateSpaceMPCController>(
                    generator, input, Q, R, states, 1, target, horizon,
                    sampleTime, limits);
            });
        result.kind = kind; result.family = "state_space";
        result.policy.evaluation = clock.mode;
        result.policy.sampleTime = clock.sampleTime;
        result.evaluation = clock.label; result.sampleTime = clock.sampleTime;
    } else if (kind == "pid") {
        const auto Kp = gainMatrix(settings.at("Kp"), 1, states, 1.0, "PID Kp");
        const auto Ki = gainMatrix(settings.at("Ki"), 1, states, 0.0, "PID Ki");
        const auto Kd = gainMatrix(settings.at("Kd"), 1, states, 0.0, "PID Kd");
        dynsyslab::core::PIDOptions pid;
        pid.derivativeFilterTime = settings.at("derivative_filter_time").number(0.0);
        pid.integratorLimit = numbers(settings.at("integrator_limit"));
        if (pid.integratorLimit.size() == 1 && states > 1) pid.integratorLimit.assign(states, pid.integratorLimit[0]);
        if (const auto cap = controlCap(settings.at("control_limit"))) {
            // Anti-windup against the actuator cap (the cap itself is applied below).
            pid.outputMinimum = cap->first; pid.outputMaximum = cap->second;
            pid.conditionalIntegration = settings.at("conditional_integration").boolean(true);
        }
        result.policy.value = controllerBank<dynsyslab::core::PIDController>(
            trajectoryCount, [Kp, Ki, Kd, states, target, pid] {
                return std::make_shared<dynsyslab::core::PIDController>(
                    Kp, Ki, Kd, states, 1, target, pid);
            });
        result.kind = kind; result.family = "state_space";
    } else if (kind == "io_delay_lqr") {
        // Echo-only ARX tracker: scalar measurement and own past commands online.
        const auto clock = resolveClock(kind, solver, settings, detail::ioDelayModelSampleTime(settings), true);
        auto io = detail::buildIoDelayControl(settings, targetSpec, system, trajectoryCount,
                                              clock.sampleTime, controlCap(settings.at("control_limit")));
        result.policy.value = std::move(io.law);
        result.policy.evaluation = clock.mode; result.policy.sampleTime = clock.sampleTime;
        result.evaluation = clock.label; result.sampleTime = clock.sampleTime;
        historyClock = true;  // the clock is bound to the model and already resolved
        result.target = std::move(io.stateTarget);
        result.observableTarget = std::move(io.measuredTarget);
        result.targetObservableNames = {io.observable};
        result.targetSpace = "observable";
        result.modelBindingEvidence = std::move(io.evidence);
        result.kind = kind; result.family = "io_delay";
    } else if (kind == "klqr" || kind == "kmpc") {
        buildKoopmanControl(kind, settings, targetSpec, solver, system, trajectoryCount, defaultSampleTime,
                            states, observableTarget, target, historyClock, result);
    } else {
        throw std::runtime_error(kind
            + " controller has not yet been ported to the native control core");
    }

    const bool sampledRequired = kind == "sdre" || kind == "pid" || (kind == "klqr"
        && (klqrIsStateful(settings) || !settings.at("bilinear_certificate").isNull()
            || !settings.at("online_measurement").isNull()));
    if (kind != "mpc" && kind != "kmpc" && !historyClock) {
        const auto clock = resolveClock(kind, solver, settings,
                                        defaultSampleTime, sampledRequired);
        result.policy.evaluation = clock.mode;
        result.policy.sampleTime = clock.sampleTime;
        result.evaluation = clock.label;
        result.sampleTime = clock.sampleTime;
    }
    if (kind == "klqr" && !settings.at("bilinear_certificate").isNull())
        detail::requireCertifiedActuation(settings, result.policy.evaluation == dynsyslab::core::ControlEvaluationMode::sampled,
                                          result.policy.sampleTime);
    // Optional actuator cap (the cap toggle), for every controller kind except MPC and
    // KMPC, whose optimizers enforce control_limit as a constraint. A scalar L clips
    // to [-L, L], a pair [low, high] clips asymmetrically. Absent, null or false: no cap.
    if (kind != "mpc" && kind != "kmpc") {
        if (const auto cap = controlCap(settings.at("control_limit"))) {
            result.policy.value = [law = std::move(result.policy.value), low = cap->first, high = cap->second](
                std::size_t trajectory, double time, const double* state) {
                return std::clamp(law(trajectory, time, state), low, high);
            };
            if (std::holds_alternative<json::Value::Object>(result.modelBindingEvidence.storage()))
                result.modelBindingEvidence.object()["control_limit"] =
                    json::Value(json::Value::Array{json::Value(cap->first), json::Value(cap->second)});
        }
    }
    if (!prescribedKind && result.policy.value)
        result.policy.value = withWindow(std::move(result.policy.value),
                                         activeWindow(settings));
    else if (!prescribedKind && !activeWindow(settings).empty())
        throw std::invalid_argument("active_time_window applies to scalar-input controllers");
    if (additive && result.policy.additiveMatrix.empty()) {  // multi-input K-ROM control set its own B
        const auto& matrix = settings.at("input_matrix").isNull() ? settings.at("B") : settings.at("input_matrix");
        result.policy.additiveMatrix = explicitMatrix(matrix, states, 1, "input_matrix");
        result.policy.inputDimension = 1;
    }
    if(liftedTarget)result.targetSpace="lifted";
    else if(!observableTarget)result.target = std::move(target);
    return result;
}

std::optional<std::pair<double, double>> controlCap(const json::Value& limit) {
    if (limit.isNull()) return std::nullopt;
    if (std::holds_alternative<bool>(limit.storage())) {
        if (!limit.boolean()) return std::nullopt;
        throw std::invalid_argument("control_limit true needs a value: give a positive number or [low, high]");
    }
    if (std::holds_alternative<json::Value::Array>(limit.storage())) {
        const auto& values = limit.array();
        if (values.size() == 1) return controlCap(values.front());
        if (values.size() != 2) throw std::invalid_argument("control_limit must be a number or [low, high]");
        const double low = values[0].number(), high = values[1].number();
        if (!std::isfinite(low) || !std::isfinite(high) || !(low < high))
            throw std::invalid_argument("control_limit [low, high] needs finite low < high");
        return std::make_pair(low, high);
    }
    const double value = limit.number(-1);
    if (!(value > 0) || !std::isfinite(value)) throw std::invalid_argument("control_limit must be positive");
    return std::make_pair(-value, value);
}

}  // namespace dynsyslab::workflow
