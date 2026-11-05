#include <dynsyslab/workflow/ControlFactory.hpp>
#include "ModelBindingOptions.hpp"
#include "HistoryControl.hpp"
#include <dynsyslab/core/HistoryObservationArchive.hpp>

#include <dynsyslab/core/ClassicalControl.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/Text.hpp>
#include <dynsyslab/workflow/ControlTargets.hpp>

#include "PrescribedControl.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

namespace dynsyslab::workflow {
namespace {

using dynsyslab::core::ScalarControlLaw;
using dynsyslab::core::StateTargetFunction;

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

std::size_t positiveCount(const json::Value& source, std::size_t fallback,
                          const char* name) {
    const double value = source.number(static_cast<double>(fallback));
    if (!std::isfinite(value) || value < 1.0 || std::floor(value) != value
        || value > 9007199254740991.0)
        throw std::invalid_argument(std::string(name) + " must be a positive integer");
    return static_cast<std::size_t>(value);
}

struct Clock {
    dynsyslab::core::ControlEvaluationMode mode{
        dynsyslab::core::ControlEvaluationMode::continuousRhs};
    std::string label{"continuous_rhs"};
    double sampleTime{};
};

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
    const std::size_t states = system.descriptor().stateDimension;
    const bool prescribedKind = kind == "waveform" || kind == "random_pulse";
    const auto& targetSpec=source.at("target");
    const auto targetSpace=targetSpec.at("space").string("state");
    if(targetSpace!="state"&&targetSpace!="observable")
        throw std::invalid_argument("control target space must be state or observable");
    const bool observableTarget=targetSpace=="observable";
    if(observableTarget&&kind!="klqr"&&kind!="kmpc")
        throw std::invalid_argument("observable-coordinate targets require Koopman control");
    StateTargetFunction target = prescribedKind || observableTarget ? StateTargetFunction{}
        : resolveControlTarget(targetSpec, states);
    const auto point = linearizationPoint(settings, fallbackState, states);
    requireScalarInput(system, point.data());

    ControllerSetup result;
    bool historyClock=false;
    if (kind == "waveform" || kind == "random_pulse") {
        auto prescribed = detail::buildPrescribedControl(kind, settings, system);
        result.kind = kind; result.family = "prescribed";
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
        result.policy.value = controllerBank<dynsyslab::core::PIDController>(
            trajectoryCount, [Kp, Ki, Kd, states, target] {
                return std::make_shared<dynsyslab::core::PIDController>(
                    Kp, Ki, Kd, states, 1, target);
            });
        result.kind = kind; result.family = "state_space";
    } else if (kind == "klqr" || kind == "kmpc") {
        const std::string path = settings.at("krom").string();
        if (path.empty())
            throw std::invalid_argument("Koopman control requires a KROM archive path");
        std::vector<std::size_t> selected;
        for (const double value : numbers(settings.at("selected_modes"))) {
            if (!std::isfinite(value) || value < 0.0 || std::floor(value) != value)
                throw std::invalid_argument("selected_modes must contain nonnegative integers");
            if (value > 9007199254740991.0)
                throw std::invalid_argument("selected_modes index is too large");
            selected.push_back(static_cast<std::size_t>(value));
        }
        auto model = std::make_shared<dynsyslab::core::KoopmanModel>(
            dynsyslab::core::loadKoopmanModel(
                dynsyslab::core::wide(path), system.configuration(), selected, modelBindingOptions(settings)));
        result.modelBindingEvidence=json::parse(model->bindingEvidenceJson());
        std::optional<dynsyslab::core::FrozenHistoryTransform> history;
        if(model->coordinates().kind==dynsyslab::core::EigenfunctionCoordinateKind::causalHistory) {
            const auto transformPath=settings.at("history_transform").string(
                model->coordinates().providerArtifact.empty()?path:model->coordinates().providerArtifact);
            history=dynsyslab::core::loadFrozenHistoryTransform(dynsyslab::core::wide(transformPath));
            const auto& identity=model->coordinates().providerIdentity;
            if(identity.rfind("sha256:",0)==0) {
                if(identity!=history->identity())throw std::invalid_argument("History transform differs from KROM training coordinates");
            } else if(settings.at("unverified_transform_reason").string().empty())
                throw std::invalid_argument("Legacy KROM has no frozen transform identity; provide unverified_transform_reason");
            if(history->retainedCoordinates().size()!=model->coordinateDimension())
                throw std::invalid_argument("History coordinate count differs from KROM");
            if(history->timeUnits()!="model-time"&&history->timeUnits()!=system.coordinates().modelTimeUnits)
                throw std::invalid_argument("History transform must use the simulator model clock; resample explicitly");
            const auto clock=resolveClock(kind,solver,settings,history->sampleInterval(),true);
            if(std::abs(clock.sampleTime-history->sampleInterval())>1e-10*history->sampleInterval())
                throw std::invalid_argument("Controller measurement clock must match the frozen history transform");
            historyClock=true;result.policy.evaluation=clock.mode;result.policy.sampleTime=clock.sampleTime;
            result.evaluation=clock.label;result.sampleTime=clock.sampleTime;
            result.modelBindingEvidence.object()["history_transform_identity"]=history->identity();
            result.modelBindingEvidence.object()["history_transform_path"]=transformPath;
            result.modelBindingEvidence.object()["unverified_transform_reason"]=settings.at("unverified_transform_reason");
        }
        if (model->inputDimension() != 1)
            throw std::invalid_argument("current native integrator supports scalar-input Koopman control");
        auto referenceSpace=dynsyslab::core::KoopmanTargetSpace::dynamicalState;
        if(observableTarget) {
            referenceSpace=dynsyslab::core::KoopmanTargetSpace::observableCoordinates;
            if(targetSpec.at("kind").string()=="from_state") {
                const auto coordinates=model->coordinateObservable();
                result.target=resolveControlTarget(targetSpec.at("state_target"),states);
                target=[stateTarget=result.target,coordinates](double time) {
                    const auto x=stateTarget(time); return coordinates.evaluateReal(x.data());
                };
            } else target=resolveControlTarget(targetSpec,model->coordinateDimension());
            result.observableTarget=target;
            result.targetObservableNames=model->coordinates().names;
            result.targetSpace="observable";
        }
        const auto Q = squareMatrix(settings.at("Q"), model->splitDimension(), "Q");
        const auto R = squareMatrix(settings.at("R"), model->inputDimension(), "R");
        const std::string representation = settings.at("control_representation").string("exact");
        const std::string mapArchive=settings.at("control_map_archive").string();
        if (!mapArchive.empty()) {
            if (representation=="exact")
                throw std::invalid_argument(
                    "an external control-map archive requires a non-exact representation name");
            model->setControlMap(representation,
                dynsyslab::core::loadApproximateControlMap(
                    dynsyslab::core::wide(mapArchive),
                    settings.at("control_map_prefix").string(),
                    *model,settings.at("unverified_control_map_reason").string()));
            result.modelBindingEvidence.object()["control_map_provenance"]=std::visit(
                [](const auto& map){return map.provenance;},model->storedControlMap(representation));
        }
        if (kind == "klqr") {
            auto controller = std::make_shared<dynsyslab::core::KoopmanLQRController>(
                model, Q, R, target, representation, referenceSpace);
            if(history)result.policy.value=detail::coordinateControllerBank(*controller,trajectoryCount);
            else result.policy.value = [controller = std::move(controller)](
                std::size_t, double time, const double* state) {
                return controller->value(time, state).front();
            };
        } else {
            const auto clock = resolveClock(kind, solver, settings, history?history->sampleInterval():defaultSampleTime, true);
            auto limits = numbers(settings.at("control_limit"));
            if (limits.empty() && !settings.at("control_limit").isNull())
                limits.push_back(settings.at("control_limit").number());
            auto controller = std::make_shared<dynsyslab::core::KoopmanMPCController>(
                model, Q, R, target,
                positiveCount(settings.at("horizon"), 20, "KMPC horizon"),
                clock.sampleTime,
                representation, std::move(limits), referenceSpace);
            if(history)result.policy.value=detail::coordinateControllerBank(*controller,trajectoryCount);
            else result.policy.value = [controller = std::move(controller)](
                std::size_t, double time, const double* state) {
                return controller->value(time, state).front();
            };
            result.policy.evaluation = clock.mode;
            result.policy.sampleTime = clock.sampleTime;
            result.evaluation = clock.label; result.sampleTime = clock.sampleTime;
        }
        if(history) {
            if(!observableTarget&&target)throw std::invalid_argument("History-coordinate control needs an observable-coordinate target");
            const auto observable=settings.at("measurement_observable").string(history->measurement());
            if(observable.empty())throw std::invalid_argument("History control requires a simulator measurement observable");
            const auto catalog=system.observableCatalog();
            if(!history->measurement().empty()&&catalog.at(observable).descriptor().name!=catalog.at(history->measurement()).descriptor().name)
                throw std::invalid_argument("Deployment measurement differs from frozen history training measurement");
            result.policy.value=detail::historyControlLaw(std::move(result.policy.value),*history,
                catalog.at(observable),trajectoryCount,settings.at("warmup_input").number(0));
        }
        result.kind = kind; result.family = "koopman";
    } else {
        throw std::runtime_error(kind
            + " controller has not yet been ported to the native control core");
    }

    const bool sampledRequired = kind == "sdre" || kind == "pid";
    if (kind != "mpc" && kind != "kmpc" && !historyClock) {
        const auto clock = resolveClock(kind, solver, settings,
                                        defaultSampleTime, sampledRequired);
        result.policy.evaluation = clock.mode;
        result.policy.sampleTime = clock.sampleTime;
        result.evaluation = clock.label;
        result.sampleTime = clock.sampleTime;
    }
    if (kind != "waveform" && kind != "random_pulse")
        result.policy.value = withWindow(std::move(result.policy.value),
                                         activeWindow(settings));
    if(!observableTarget)result.target = std::move(target);
    return result;
}

}  // namespace dynsyslab::workflow
