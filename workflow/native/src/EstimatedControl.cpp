#include "KlqrAugmentation.hpp"
#include "EstimatedControl.hpp"
#include "ModelBindingOptions.hpp"
#include "BalancedTruncationSettings.hpp"
#include "BilinearCertificateBinding.hpp"
#include <dynsyslab/workflow/KoopmanFeedback.hpp>
#include <dynsyslab/workflow/EstimatorFactory.hpp>
#include <dynsyslab/workflow/ControlTargets.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/ScientificCompatibility.hpp>
#include <dynsyslab/core/Text.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>

namespace dynsyslab::workflow::detail {
namespace {
std::vector<double> numbers(const json::Value& v) {
    std::vector<double> result;
    for(const auto& item:v.array())result.push_back(item.number());
    return result;
}
std::vector<double> square(const json::Value& v,std::size_t n,const char* label) {
    auto values=numbers(v);
    if(values.size()==n*n)return values;
    std::vector<double> result(n*n,0);
    if(!values.empty()&&values.size()!=n)
        throw std::invalid_argument(std::string(label)+" must be scalar, diagonal or square");
    for(std::size_t k=0;k<n;++k)result[k*n+k]=values.empty()?v.number(1):values[k];
    return result;
}
struct Trace {
    std::vector<double> time,lifted,output,measured,predicted,innovation,requested,applied,previous;
    void append(const KoopmanFeedbackSample& sample,double actual) {
        time.push_back(sample.estimate.time);
        const auto add=[](auto& a,const auto& b){a.insert(a.end(),b.begin(),b.end());};
        add(lifted,sample.estimate.lifted);add(output,sample.estimate.output);
        add(measured,sample.measured);add(predicted,sample.estimate.measurement);
        add(innovation,sample.estimate.innovation);add(requested,sample.commandedInput);
        add(previous,sample.previousAppliedInput);applied.push_back(actual);
    }
    void store(core::NpzRecord& archive) const {
        const auto rows=time.size();
        archive.arrays["estimator_time"]=core::NpyArray::fromDoubles({rows},time);
        const auto save=[&](const char* name,const std::vector<double>& values) {
            if(!values.empty())archive.arrays[name]=core::NpyArray::fromDoubles(
                {rows,values.size()/rows},values);
        };
        save("lifted_estimate",lifted);save("state_estimate",output);
        save("measurement_values",measured);save("measurement_estimate",predicted);
        save("innovation",innovation);save("commanded_input",requested);
        save("applied_input",applied);save("previous_applied_input",previous);
    }
};
struct Stream {
    std::unique_ptr<KoopmanFeedback> feedback;
    double applied{};
    Trace trace;
    SensorDifferentiator derivative;
};
}

ControllerSetup buildEstimatedController(const json::Value& source,
    const json::Value& solver,const core::ODESystem& system,
    std::size_t count,double defaultStep,
    std::function<double(std::size_t,double,const double*)> appliedOverride) {
    if(source.at("kind").string()!="klqr"||!count)
        throw std::invalid_argument("Online KLQE composition currently requires KLQR");
    const auto& settings=source.at("settings");
    const auto& estimator=settings.at("estimator");
    if(estimator.at("kind").string("klqe")!="klqe")
        throw std::invalid_argument("Unknown online Koopman estimator");
    const auto evaluation=solver.at("control_evaluation").string("automatic");
    if(evaluation!="automatic"&&evaluation!="sampled")
        throw std::invalid_argument("KLQE feedback requires sampled control evaluation");
    const double dt=solver.at("control_sample_time").number(
        settings.at("sample_time").number(defaultStep));
    if(!std::isfinite(dt)||dt<=0)
        throw std::invalid_argument("KLQE feedback requires a positive sample time");
    const auto path=settings.at("krom").string();
    if(path.empty())throw std::invalid_argument("KLQE feedback requires a KROM archive");
    std::vector<std::size_t> selected;
    for(double v:numbers(settings.at("selected_modes"))) {
        if(!std::isfinite(v)||v<0||v!=std::floor(v)||v>9007199254740991.0)
            throw std::invalid_argument("selected_modes must contain nonnegative integers");
        selected.push_back(static_cast<std::size_t>(v));
    }
    auto model=core::loadKoopmanModel(core::wide(path),system.configuration(),
        selected,modelBindingOptions(settings));
    if(model.inputDimension()!=1)
        throw std::invalid_argument("ODE solver adapter currently requires one input");
    const auto representation=settings.at("control_representation").string();
    if(representation.empty()||representation=="exact")
        throw std::invalid_argument("KLQE feedback requires an explicit fitted control_representation");
    const auto mapPath=settings.at("control_map_archive").string();
    if(!mapPath.empty())model.setControlMap(representation,core::loadApproximateControlMap(
        core::wide(mapPath),settings.at("control_map_prefix").string(),model,
        settings.at("unverified_control_map_reason").string()));
    // KLQE binding shared with the Estimator module (design, sensors, initial estimate, intercepts).
    auto bound=bindKlqe(estimator,model,system);
    bound.sampleTime=dt;
    const auto& measurementName=bound.measurement;
    const auto& outputName=bound.output;
    const auto n=model.splitDimension(),p=model.reconstruction(measurementName).outputs;
    const auto design=bound.design;
    auto sensor=std::move(bound.sensor);
    const auto initial=bound.initial;
    const auto offsetY=bound.measurementOffset;
    const auto offsetX=bound.outputOffset;
    const auto feedthrough=bound.measurementFeedthrough;
    const int feedthroughScale=bound.feedthroughScaleOutput;
    const bool lift=bound.liftSensors;
    const auto liftModel=std::make_shared<const core::KoopmanModel>(model);
    const auto Q=square(settings.at("Q"),n,"Q"),R=square(settings.at("R"),1,"R");
    const auto cap=controlCap(settings.at("control_limit"));
    const double lowLimit=cap?cap->first:-(std::numeric_limits<double>::infinity)();
    const double highLimit=cap?cap->second:(std::numeric_limits<double>::infinity)();
    const auto window=numbers(settings.at("active_time_window"));
    if(!window.empty()&&(window.size()!=2||!std::isfinite(window[0])
        ||!std::isfinite(window[1])||window[0]>window[1]))
        throw std::invalid_argument("active_time_window must be [start, stop]");

    ControllerSetup result;
    core::StateTargetFunction target;
    auto targetSpace=core::KoopmanTargetSpace::realLiftedValues;
    const auto& targetSpec=source.at("target");
    result.targetSpace=targetSpec.at("space").string("lifted");
    if(!targetSpec.isNull()) {
        std::size_t dimension=n;
        if(result.targetSpace=="state") {
            targetSpace=core::KoopmanTargetSpace::dynamicalState;
            dimension=model.stateDimension();
            if(model.requiresExternalCoordinates())
                throw std::invalid_argument("History KROM feedback requires lifted or observable targets");
        } else if(result.targetSpace=="observable") {
            targetSpace=core::KoopmanTargetSpace::observableCoordinates;
            dimension=model.coordinateDimension();
        } else if(result.targetSpace!="lifted")
            throw std::invalid_argument("Estimated KLQR target space must be state, observable or lifted");
        target=resolveControlTarget(targetSpec,dimension);
        if(result.targetSpace=="state")result.target=target;
        else {
            result.observableTarget=target;
            if(result.targetSpace=="observable")result.targetObservableNames=model.coordinates().names;
            else for(std::size_t k=0;k<n;++k)
                result.targetObservableNames.push_back("z"+std::to_string(k));
        }
    }
    const auto controllerReduction=parseBalancedTruncation(settings.at("balanced_truncation"));
    std::optional<BoundBilinearCertificate> certificate;
    if(!settings.at("bilinear_certificate").isNull()) {
        certificate=bindBilinearCertificate(settings,model,representation);
        requireCertifiedActuation(settings,true,dt);
    }
    auto streams=std::make_shared<std::vector<Stream>>(count);
    for(auto& stream:*streams) {
        stream.feedback=std::make_unique<KoopmanFeedback>(
            model,design,representation,initial,Q,R,target,targetSpace,offsetY,offsetX);
        stream.feedback->setMeasurementFeedthrough(feedthrough,feedthroughScale);
        if(controllerReduction)stream.feedback->setBalancedReduction(*controllerReduction);
        // Integral action / PID on a KROM readout of the estimate z_hat (e.g. r_hat).
        applyKlqrAugmentations(stream.feedback->controller(),model,settings);
        if(certificate)stream.feedback->controller().setCertifiedFeedback(certificate->feedback);
        stream.derivative=SensorDifferentiator(bound,std::size_t(&stream-&streams->front()));
    }
    const bool shadow=static_cast<bool>(appliedOverride);
    result.policy.value=[streams,sensor=std::move(sensor),lowLimit,highLimit,window,lift,liftModel,
        appliedOverride=std::move(appliedOverride)](
        std::size_t trajectory,double time,const double* state) {
        auto& stream=streams->at(trajectory);
        // This is the only truth-state access: it simulates the sensor.
        auto y=sensor->evaluate(state,stream.applied);
        stream.derivative.append(y,time);  // numerically differentiated sensor outputs
        if(lift)y=liftModel->evaluateCoordinatesSplit(y.data());  // psi of the sensor values
        const auto sample=stream.feedback->sample(time,y,{stream.applied});
        double applied=std::clamp(sample.commandedInput.front(),lowLimit,highLimit);
        // Gating happens AFTER observation, so estimator memory is never frozen
        // outside the control window and the actual zero input is accounted for.
        if(!window.empty()&&(time<window[0]||time>window[1]))applied=0;
        // applied_controller (shadow mode): the plant receives the reference drive instead; the
        // KLQE predicts with that drive and the KLQR command is only archived (commanded_input).
        if(appliedOverride)applied=appliedOverride(trajectory,time,state);
        stream.trace.append(sample,applied);stream.applied=applied;
        return applied;
    };
    result.appendDiagnostics=[streams,design,n,p](std::size_t trajectory,core::NpzRecord& archive) {
        streams->at(trajectory).trace.store(archive);
        storeKlqrReduction(streams->at(trajectory).feedback->controller(),archive);
        archive.arrays["estimator_gain"]=core::NpyArray::fromDoubles({n,p},design.solution.gain);
        archive.arrays["observer_poles"]=core::NpyArray::fromComplex({n},design.solution.observerPoles);
        archive.arrays["measurement_matrix"]=core::NpyArray::fromDoubles({p,n},design.measurementMatrix);
        if(design.outputDimension)archive.arrays["state_reconstruction_matrix"]
            =core::NpyArray::fromDoubles({design.outputDimension,n},design.stateReconstructionMatrix);
    };
    result.kind="klqr";result.family="koopman";result.evaluation="sampled";
    result.sampleTime=dt;result.policy.sampleTime=dt;
    result.policy.evaluation=core::ControlEvaluationMode::sampled;
    result.modelBindingEvidence=json::parse(model.bindingEvidenceJson());
    auto& evidence=result.modelBindingEvidence.object();
    evidence["feedback_architecture"]=shadow?"klqe-klqr-shadow-of-applied-controller/v1"
                                            :"klqe-estimated-lifted-klqr/v1";
    if(shadow)evidence["applied_controller"]=settings.at("applied_controller");
    evidence["model_identity"]=core::koopmanModelIdentity(model);
    evidence["measurement_reconstruction"]=measurementName;
    evidence["output_reconstruction"]=outputName;
    evidence["observer_integration"]="exact-zoh-previous-measurement-and-applied-input";
    evidence["control_map_provenance"]=std::visit(
        [](const auto& map){return map.provenance;},model.storedControlMap(representation));
    evidence["estimator_clock"]="estimator_time";
    evidence["sensor_feedthrough"]=estimator.at("sensor_feedthrough");
    evidence["sensor_input_timing"]="previous-applied-input-before-new-command";
    if(klqrIsStateful(settings))evidence["klqr_augmentations"]=klqrAugmentationEvidence(settings);
    if(certificate) {
        auto certified=certificate->evidence;
        certified.object()["output_feedback"]="not certified: the design assumes z is known, and here u uses the KLQE estimate";
        certified.object()["klqe_prediction"]=certified.at("operators").string()=="fit"
            ?"KLQE propagates the stored continuous bilinear map; the gain was designed on B0 and B_i refitted from held-input records, so the two input models differ"
            :"KLQE propagates the same frozen generator and stored bilinear map, exact ZOH for each held input";
        evidence["bilinear_certificate"]=std::move(certified);
    }
    if(!bound.reduction.isNull())evidence["estimator_balanced_truncation"]=bound.reduction;
    if(!bound.processCovariance.isNull())evidence["estimator_process_covariance"]=bound.processCovariance;
    if(controllerReduction)evidence["controller_balanced_truncation"]=balancedTruncationRequest(*controllerReduction);
    return result;
}
} // namespace dynsyslab::workflow::detail
