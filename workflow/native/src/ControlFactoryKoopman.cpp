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
#include "KromDomainPolicy.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>


namespace dynsyslab::workflow::control_factory_detail {

void buildKoopmanControl(const std::string& kind, const json::Value& settings, const json::Value& targetSpec,
                         const json::Value& solver, const dynsyslab::core::ODESystem& system,
                         std::size_t trajectoryCount, double defaultSampleTime, std::size_t states,
                         bool observableTarget, dynsyslab::core::StateTargetFunction& target,
                         bool& historyClock, ControllerSetup& result) {
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
    if(model->coordinates().computableFromMeasurements()) {
        auto map=postprocess_detail::loadMeasurementMap(*model,dynsyslab::core::wide(path),settings);
        const auto transformPath=map.path;
        history=std::move(map.transform);
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
    }
    // actuation "additive_state" with an explicit n x m B (settings.B or input_matrix): the exact
    // map becomes grad psi(x) B and KLQR/KMPC drive m inputs through the plant's xdot = f(x) + B u.
    if(settings.at("actuation").string("system_input")=="additive_state") {
        const auto& matrix=settings.at("input_matrix").isNull()?settings.at("B"):settings.at("input_matrix");
        model->setAdditiveInputMatrix(numbers(matrix));
        if(model->additiveInputMatrix().empty())
            throw std::invalid_argument("additive_state Koopman control needs B (state dimension x inputs)");
    }
    const std::size_t inputs=model->inputDimension();
    if (inputs != 1 && (history || !settings.at("online_measurement").isNull()))
        throw std::invalid_argument("multi-input Koopman control needs state feedback: "
                                    "no causal-history or online-measurement coordinates");
    // One law for both input shapes: a scalar law, or a vector law with additive actuation.
    const auto assign=[&](auto evaluate) {
        if(inputs==1) {
            result.policy.value=[evaluate](std::size_t k,double time,const double* state) {
                return evaluate(k,time,state).front();
            };
            return;
        }
        result.policy.vector=[evaluate](std::size_t k,double time,const double* state,double* input) {
            const auto u=evaluate(k,time,state); std::copy(u.begin(),u.end(),input);
        };
        result.policy.inputDimension=inputs;
        result.policy.additiveMatrix=model->additiveInputMatrix();
    };
    auto referenceSpace=dynsyslab::core::KoopmanTargetSpace::dynamicalState;
    const auto historyMeasurement=[&]() {
        const auto name=settings.at("measurement_observable").string(history->measurement());
        if(name.empty())throw std::invalid_argument("History control requires a simulator measurement observable");
        return system.observableCatalog().at(name);
    };
    if(observableTarget) {
        referenceSpace=dynsyslab::core::KoopmanTargetSpace::observableCoordinates;
        if(targetSpec.at("kind").string()=="from_state") {
            result.target=resolveControlTarget(targetSpec.at("state_target"),states);
            if(history) {
                // xi_ref: frozen projection of the measured observable along the state target.
                target=detail::historyCoordinateTarget(result.target,*history,historyMeasurement());
            } else {
                const auto coordinates=model->coordinateObservable();
                target=[stateTarget=result.target,coordinates](double time) {
                    const auto x=stateTarget(time); return coordinates.evaluateReal(x.data());
                };
            }
        } else target=resolveControlTarget(targetSpec,model->coordinateDimension());
        result.observableTarget=target;
        result.targetObservableNames=model->coordinates().names;
        result.targetSpace="observable";
    } else if(targetSpec.at("space").string()=="lifted") {
        // A reference z_ref(t) given directly in the real split lifted coordinates.
        referenceSpace=dynsyslab::core::KoopmanTargetSpace::realLiftedValues;
        target=resolveControlTarget(targetSpec,model->splitDimension());
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
                *model, settings.at("control_map_origin").string() == "external"));
        result.modelBindingEvidence.object()["control_map_provenance"]=std::visit(
            [](const auto& map){return map.provenance;},model->storedControlMap(representation));
    }
    const auto reduction=parseBalancedTruncation(settings.at("balanced_truncation"));
    if(reduction&&kind!="klqr")
        throw std::invalid_argument("balanced_truncation applies to klqr (and klqe), not "+kind);
    if(reduction)result.modelBindingEvidence.object()["controller_balanced_truncation"]=
        balancedTruncationRequest(*reduction);
    const bool historyExact=representation=="history_exact";
    std::optional<detail::BoundBilinearCertificate> certificate;
    if(!settings.at("bilinear_certificate").isNull()) {
        if(kind!="klqr"||historyExact)
            throw std::invalid_argument("bilinear_certificate applies to klqr with a stored bilinear map");
        certificate=detail::bindBilinearCertificate(settings,*model,representation);
        result.modelBindingEvidence.object()["bilinear_certificate"]=certificate->evidence;
    }
    if(historyExact&&(!history||kind!="klqr"))
        throw std::invalid_argument("control_representation history_exact needs a causal-history KROM and klqr");
    if (kind == "klqr" && historyExact) {
        auto controller = std::make_shared<dynsyslab::core::KoopmanLQRController>(
            model, Q, R, target, "external", referenceSpace);
        if(reduction)controller->setBalancedReduction(*reduction);
        std::shared_ptr<const dynsyslab::core::ODESystem> plant=
            dynsyslab::core::createSystem(system.configuration());
        result.policy.value=detail::historyExactControlLaw(std::move(controller),model,*history,
            historyMeasurement(),std::move(plant),trajectoryCount,settings.at("warmup_input").number(0));
        result.modelBindingEvidence.object()["control_map_provenance"]=
            "history_exact: grad_xi psi * newest-sample weight * L_g(measurement), from simulator state";
    } else if (kind == "klqr") {
        auto controller = std::make_shared<dynsyslab::core::KoopmanLQRController>(
            model, Q, R, target, representation, referenceSpace);
        if(reduction)controller->setBalancedReduction(*reduction);
        applyKlqrAugmentations(*controller, *model, settings);
        if(certificate)controller->setCertifiedFeedback(certificate->feedback);
        if(controller->stateful())
            result.modelBindingEvidence.object()["klqr_augmentations"]=klqrAugmentationEvidence(settings);
        if(const auto& online=settings.at("online_measurement");!online.isNull()) {
            // Signal-only feedback: coordinates (y, y') from one sensor and a causal differentiator.
            if(history||controller->stateful()||representation=="exact"||model->coordinateDimension()!=2)
                throw std::invalid_argument("online_measurement needs a two-coordinate KROM, a stored input map and a stateless KLQR");
            const auto clock=resolveClock(kind,solver,settings,defaultSampleTime,true);
            dynsyslab::core::DerivativeSpec derivative;derivative.method="backward2";
            if(!online.at("derivative").isNull())derivative=parseDerivativeSpec(online.at("derivative"));
            result.policy.value=detail::measuredJetControlLaw(*controller,
                system.observableCatalog().at(online.at("observable").string("far_field_acoustic_echo")),
                derivative,clock.sampleTime,trajectoryCount,
                static_cast<std::size_t>(online.at("warmup_samples").number(3)),settings.at("warmup_input").number(0));
            result.modelBindingEvidence.object()["online_measurement"]=online;
        } else if(history)result.policy.value=detail::coordinateControllerBank(*controller,trajectoryCount);
        else if(controller->stateful()) {
            // Integral/PID states: one controller copy per trajectory, sampled clock (below).
            auto bank=std::make_shared<std::vector<dynsyslab::core::KoopmanLQRController>>(
                trajectoryCount,*controller);
            assign([bank](std::size_t trajectory,double time,const double* state) {
                return bank->at(trajectory).value(time,state);
            });
        } else {
            if(reduction)result.appendDiagnostics=[controller](std::size_t,dynsyslab::core::NpzRecord& archive) {
                storeKlqrReduction(*controller,archive);
            };
            assign([controller](std::size_t, double time, const double* state) {
                return controller->value(time, state);
            });
        }
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
        else assign([controller](std::size_t, double time, const double* state) {
            return controller->value(time, state);
        });
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
        if(!historyExact)  // history_exact already built its own measurement law
            result.policy.value=detail::historyControlLaw(std::move(result.policy.value),*history,
                catalog.at(observable),trajectoryCount,settings.at("warmup_input").number(0));
    }
    // target_domain_policy: the plant state and the state target checked against the K-ROM's
    // training domain on every evaluation (error | warn | ignore; KromDomainPolicy.cpp).
    detail::applyTargetDomainPolicy(result, settings, path, states,
        targetSpec.at("space").string("state")=="state"?target:dynsyslab::core::StateTargetFunction{},
        trajectoryCount);
    result.kind = kind; result.family = "koopman";
}

}  // namespace dynsyslab::workflow::control_factory_detail
