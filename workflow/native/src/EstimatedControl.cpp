#include "EstimatedControl.hpp"
#include "ModelBindingOptions.hpp"
#include <dynsyslab/workflow/KoopmanFeedback.hpp>
#include <dynsyslab/workflow/ControlTargets.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/ScientificCompatibility.hpp>
#include <dynsyslab/core/Text.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
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
};
}

ControllerSetup buildEstimatedController(const json::Value& source,
    const json::Value& solver,const core::ODESystem& system,
    std::size_t count,double defaultStep) {
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
    const auto measurementName=estimator.at("measurement").string();
    const auto outputName=estimator.at("output").string();
    if(measurementName.empty())
        throw std::invalid_argument("KLQE requires a named measurement reconstruction");
    const auto n=model.splitDimension(),p=model.reconstruction(measurementName).outputs;
    const auto design=core::designKoopmanLQEFromMeasurements(model,measurementName,
        square(estimator.at("process_covariance"),n,"process covariance"),
        square(estimator.at("measurement_covariance"),p,"measurement covariance"),outputName);
    std::vector<core::Observable> sensors;
    const auto catalog=system.observableCatalog();
    if(!estimator.at("measurement_observables").array().empty()) {
        for(const auto& name:estimator.at("measurement_observables").array())
            sensors.push_back(catalog.at(name.string()));
    } else sensors.push_back(catalog.at(
        estimator.at("measurement_observable").string(measurementName)));
    std::size_t sensorWidth=0;
    for(const auto& sensor:sensors)sensorWidth+=sensor.outputDimension();
    if(sensorWidth!=p)
        throw std::invalid_argument("Simulated sensor outputs differ from C_y row count");
    auto initial=numbers(estimator.at("initial_lifted"));
    if(initial.empty())initial.assign(n,0);
    const auto offsetY=numbers(estimator.at("measurement_offset"));
    const auto offsetX=numbers(estimator.at("output_offset"));
    const auto Q=square(settings.at("Q"),n,"Q"),R=square(settings.at("R"),1,"R");
    double limit=(std::numeric_limits<double>::infinity)();
    if(!settings.at("control_limit").isNull()) {
        limit=settings.at("control_limit").number(-1);
        if(!std::isfinite(limit)||limit<=0)
            throw std::invalid_argument("KLQE scalar control_limit must be positive");
    }
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
    auto streams=std::make_shared<std::vector<Stream>>(count);
    for(auto& stream:*streams)stream.feedback=std::make_unique<KoopmanFeedback>(
        model,design,representation,initial,Q,R,target,targetSpace,offsetY,offsetX);
    result.policy.value=[streams,sensors=std::move(sensors),limit,window](
        std::size_t trajectory,double time,const double* state) {
        auto& stream=streams->at(trajectory);
        std::vector<double> y;
        // This is the only truth-state access: it simulates the sensor.
        for(const auto& sensor:sensors) {
            const auto values=sensor.evaluateReal(state);y.insert(y.end(),values.begin(),values.end());
        }
        const auto sample=stream.feedback->sample(time,y,{stream.applied});
        double applied=std::clamp(sample.commandedInput.front(),-limit,limit);
        // Gating happens AFTER observation, so estimator memory is never frozen
        // outside the control window and the actual zero input is accounted for.
        if(!window.empty()&&(time<window[0]||time>window[1]))applied=0;
        stream.trace.append(sample,applied);stream.applied=applied;
        return applied;
    };
    result.appendDiagnostics=[streams,design,n,p](std::size_t trajectory,core::NpzRecord& archive) {
        streams->at(trajectory).trace.store(archive);
        archive.arrays["estimator_gain"]=core::NpyArray::fromDoubles({n,p},design.solution.gain);
        archive.arrays["measurement_matrix"]=core::NpyArray::fromDoubles({p,n},design.measurementMatrix);
        if(design.outputDimension)archive.arrays["state_reconstruction_matrix"]
            =core::NpyArray::fromDoubles({design.outputDimension,n},design.stateReconstructionMatrix);
    };
    result.kind="klqr";result.family="koopman";result.evaluation="sampled";
    result.sampleTime=dt;result.policy.sampleTime=dt;
    result.policy.evaluation=core::ControlEvaluationMode::sampled;
    result.modelBindingEvidence=json::parse(model.bindingEvidenceJson());
    auto& evidence=result.modelBindingEvidence.object();
    evidence["feedback_architecture"]="klqe-estimated-lifted-klqr/v1";
    evidence["model_identity"]=core::koopmanModelIdentity(model);
    evidence["measurement_reconstruction"]=measurementName;
    evidence["output_reconstruction"]=outputName;
    evidence["observer_integration"]="exact-zoh-previous-measurement-and-applied-input";
    evidence["control_map_provenance"]=std::visit(
        [](const auto& map){return map.provenance;},model.storedControlMap(representation));
    evidence["estimator_clock"]="estimator_time";
    return result;
}
} // namespace dynsyslab::workflow::detail
