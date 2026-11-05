#include <dynsyslab/workflow/EstimationPostprocess.hpp>
#include "ObservationSource.hpp"
#include "ModelBindingOptions.hpp"

#include <dynsyslab/core/Text.hpp>

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Estimation.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/System.hpp>

#include <algorithm>
#include <chrono>
#include <cctype>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace dynsyslab::workflow {
namespace {

std::string escaped(const std::string& value) {
    std::string result;
    for (const char item:value) {
        if (item=='"' || item=='\\') result.push_back('\\');
        if (item=='\n') result+="\\n"; else if (item!='\r') result.push_back(item);
    }
    return result;
}


std::vector<double> diagonal(std::size_t size,double scalar,const json::Value& values) {
    std::vector<double> result(size*size,0.0);
    if (!values.array().empty()) {
        if (values.array().size()!=size) throw std::invalid_argument("covariance diagonal has wrong size");
        for (std::size_t i=0;i<size;++i) result[i*size+i]=values.array()[i].number();
    } else for (std::size_t i=0;i<size;++i) result[i*size+i]=scalar;
    return result;
}

std::string complexJson(const std::vector<std::complex<double>>& values) {
    std::ostringstream out; out<<'[';
    for (std::size_t i=0;i<values.size();++i) {
        if (i) out<<',';
        out<<"{\"real\":"<<std::setprecision(17)<<values[i].real()
           <<",\"imag\":"<<values[i].imag()<<'}';
    }
    out<<']'; return out.str();
}


std::vector<std::size_t> selectedModes(const json::Value& settings) {
    std::vector<std::size_t> result;
    for (const auto& value:settings.at("selected_modes").array())
        result.push_back(static_cast<std::size_t>(value.number()));
    return result;
}

dynsyslab::core::ObservabilityOptions observabilityOptions(
    const json::Value& settings,std::size_t stateDimension) {
    dynsyslab::core::ObservabilityOptions result;
    const std::string method=settings.at("observability_method").string("matrix");
    if (method=="matrix")
        result.method=dynsyslab::core::ObservabilityMethod::classicalMatrix;
    else if (method=="discrete-gramian")
        result.method=dynsyslab::core::ObservabilityMethod::discreteGramian;
    else if (method=="continuous-gramian")
        result.method=dynsyslab::core::ObservabilityMethod::continuousFiniteGramian;
    else if (method=="continuous-infinite-gramian")
        result.method=dynsyslab::core::ObservabilityMethod::continuousInfiniteGramian;
    else throw std::invalid_argument("unknown observability method: "+method);
    const double blocks=settings.at("observability_blocks").number();
    if (blocks!=0.0) {
        if (!std::isfinite(blocks) || blocks<1.0 || std::floor(blocks)!=blocks
            || blocks>9007199254740991.0)
            throw std::invalid_argument("observability_blocks must be a positive integer");
        result.blockCount=static_cast<std::size_t>(blocks);
    } else result.blockCount=stateDimension;
    result.horizon=settings.at("observability_horizon").number(1.0);
    result.sampleTime=settings.at("observability_sample_time").number();
    result.relativeTolerance=settings.at("observability_tolerance").number(1.0e-10);
    result.maximumAbsoluteEntry=
        settings.at("observability_maximum_entry").number(1.0e150);
    return result;
}

std::string numberJson(double value) {
    if (!std::isfinite(value)) return "null";
    std::ostringstream out; out<<std::setprecision(17)<<value; return out.str();
}

}  // namespace

std::filesystem::path runEstimationPostprocess(
    const std::string& method,const std::filesystem::path& sourcePath,
    const json::Value& settings,const std::filesystem::path& root,const std::string& label) {
    const auto source=dynsyslab::core::readNpz(sourcePath);
    const std::string measurement=settings.at("measurement").string("echo");
    const auto arrayName=settings.at("measurement_array").string(
        source.contains("measurement_values")?"measurement_values":measurement);
    const auto observed=postprocess_detail::loadObservationSource(sourcePath,arrayName,measurement);
    const auto config=postprocess_detail::observationSystem(observed);
    const std::string kromPath=settings.at("krom").string();
    if (kromPath.empty()) throw std::invalid_argument("estimation requires a KROM path");
    auto model=config?dynsyslab::core::loadKoopmanModel(
            dynsyslab::core::wide(kromPath),*config,selectedModes(settings),modelBindingOptions(settings)):
        dynsyslab::core::loadKoopmanModel(dynsyslab::core::wide(kromPath),selectedModes(settings));
    if(!config&&model.hasPhysicalSystem()&&settings.at("unverified_legacy_reason").string().empty())
        throw std::invalid_argument("Measurement source lacks physical provenance; provide unverified_legacy_reason for a physical KROM");
    const std::string stateName=settings.at("output_reconstruction").string(
        settings.at("state_reconstruction").string(model.hasReconstruction("state")?"state":""));
    const std::size_t n=model.splitDimension();
    if (!model.hasReconstruction(measurement) || (!stateName.empty()&&!model.hasReconstruction(stateName)))
        throw std::invalid_argument("requested KROM reconstruction map is unavailable");
    const std::size_t p=model.reconstruction(measurement).outputs;
    if(observed.dimension!=p)throw std::invalid_argument("Measurement components differ from the KROM reconstruction");
    const auto Q=diagonal(n,settings.at("process_variance").number(1.0),
                          settings.at("process_diagonal"));
    const auto R=diagonal(p,settings.at("measurement_variance").number(1.0e-4),
                          settings.at("measurement_diagonal"));
    const auto design=dynsyslab::core::designKoopmanLQEFromMeasurements(
        model,measurement,Q,R,stateName,settings.at("require_detectable").boolean(true));
    const auto observability=dynsyslab::core::analyzeObservability(
        design.generator,design.measurementMatrix,n,p,
        observabilityOptions(settings,n));
    dynsyslab::core::NpzRecord output;
    if (!observability.matrix.empty())
        output.arrays["observability_matrix"]=dynsyslab::core::NpyArray::fromDoubles(
            {observability.outputDimension*observability.blockCount,n},observability.matrix);
    output.arrays["observability_gramian"]=dynsyslab::core::NpyArray::fromDoubles(
        {n,n},observability.gramian);
    output.arrays["observability_singular_values"]=dynsyslab::core::NpyArray::fromDoubles(
        {observability.singularValues.size()},observability.singularValues);
    output.arrays["observable_directions"]=dynsyslab::core::NpyArray::fromDoubles(
        {n,n},observability.rightDirections);
    output.arrays["estimator_gain"]=dynsyslab::core::NpyArray::fromDoubles(
        {n,p},design.solution.gain);
    output.arrays["estimator_covariance"]=dynsyslab::core::NpyArray::fromDoubles(
        {n,n},design.solution.covariance);
    output.arrays["observer_poles"]=dynsyslab::core::NpyArray::fromComplex(
        {n},design.solution.observerPoles);
    output.arrays["measurement_matrix"]=dynsyslab::core::NpyArray::fromDoubles(
        {p,n},design.measurementMatrix);
    if(design.outputDimension)output.arrays["state_reconstruction_matrix"]=dynsyslab::core::NpyArray::fromDoubles(
        {design.outputDimension,n},design.stateReconstructionMatrix);
    std::string details="\"source_krom\":\""+escaped(kromPath)
        +"\",\"measurement_observable\":\""+escaped(measurement)
        +"\",\"state_reconstruction\":\""+escaped(stateName)
        +"\",\"observability_method\":\""
        +dynsyslab::core::observabilityMethodName(observability.method)
        +"\",\"observability_rank\":"+std::to_string(observability.rank)
        +",\"observability_horizon\":"+numberJson(observability.horizon)
        +",\"observability_condition_number\":"
        +numberJson(observability.conditionNumber)
        +",\"observability_effective_rank\":"
        +numberJson(observability.effectiveRank)
        +",\"observability_rank_tolerance\":"
        +numberJson(observability.rankTolerance)
        +",\"lifted_dimension\":"+std::to_string(n)
        +",\"detectable\":"+(design.solution.detectable?"true":"false")
        +",\"observer_poles\":"+complexJson(design.solution.observerPoles)
        +",\"model_binding_evidence\":"+model.bindingEvidenceJson()
        +",\"unverified_measurement_reason\":"+json::dump(settings.at("unverified_legacy_reason"));
    if (method=="klqe") {
        const auto sourceMeta=source.metadataJson().empty()?json::Value{}:json::parse(source.metadataJson());
        const auto space=sourceMeta.at("value_space").string();
        const bool stateSpace=space.empty()||space=="physical-state"||space=="dynamical-state";
        const auto stateKey=source.contains("states")?"states":"values";
        const bool hasTruth=config&&stateSpace&&model.physicalStateDimension()&&source.contains(stateKey)
            &&source.at(stateKey).shape==std::vector<std::size_t>{observed.samples,*model.physicalStateDimension()};
        std::vector<double> initial,truth;
        for(const auto& value:settings.at("initial_lifted_values").array())initial.push_back(value.number());
        if(initial.empty()&&!settings.at("initial_state").isNull()) {
            std::vector<double> state;for(const auto& value:settings.at("initial_state").array())state.push_back(value.number());
            if(state.size()!=model.stateDimension())throw std::invalid_argument("Initial state dimension differs from KROM");
            initial=model.evaluateSplit(state.data());
        }
        if(initial.empty()&&hasTruth&&!model.requiresExternalCoordinates()) {
            const auto system=dynsyslab::core::createSystem(*config);
            initial=model.evaluateSplit(system->initialState().data());
        }
        if(initial.empty())initial.assign(n,0.0);
        if(stateName=="state"&&hasTruth) {
            const auto physical=postprocess_detail::loadSource(sourcePath);
            if(physical.time!=observed.time)throw std::invalid_argument("Truth and measurement clocks differ; explicitly normalize measurements first");
            truth=physical.values;
        }
        const auto& measurements=observed.values;
        const auto run=dynsyslab::core::runKoopmanLQEFromLifted(
            model,design,observed.time,measurements,initial,truth);
        output.arrays["time"]=dynsyslab::core::NpyArray::fromDoubles({run.samples},run.time);
        output.arrays["values"]=dynsyslab::core::NpyArray::fromDoubles(
            {run.samples,run.stateDimension?run.stateDimension:run.liftedDimension},
            run.stateDimension?run.stateEstimate:run.liftedEstimate);
        if(!truth.empty())output.arrays["true_states"]=dynsyslab::core::NpyArray::fromDoubles(
            {run.samples,run.stateDimension},truth);
        output.arrays["lifted_estimate"]=dynsyslab::core::NpyArray::fromDoubles(
            {run.samples,run.liftedDimension},run.liftedEstimate);
        output.arrays["measurement_values"]=dynsyslab::core::NpyArray::fromDoubles(
            {run.samples,run.measurementDimension},measurements);
        output.arrays["measurement_estimate"]=dynsyslab::core::NpyArray::fromDoubles(
            {run.samples,run.measurementDimension},run.measurementEstimate);
        output.arrays["innovation"]=dynsyslab::core::NpyArray::fromDoubles(
            {run.samples,run.measurementDimension},run.innovation);
        if(!run.stateError.empty())output.arrays["state_error"]=dynsyslab::core::NpyArray::fromDoubles(
            {run.samples,run.stateDimension},run.stateError);
        output.arrays["state_rmse"]=dynsyslab::core::NpyArray::fromDoubles(
            {run.stateDimension},run.stateRmse);
        output.arrays["innovation_rmse"]=dynsyslab::core::NpyArray::fromDoubles(
            {run.measurementDimension},run.innovationRmse);
        details+=",\"integration\":\"exact-zero-order-held-measurement\"";
        details+=",\"value_space\":\""+std::string(run.stateDimension?"reconstructed-output":"real-lifted")+"\"";
    }
    output.arrays["metadata_json"]=dynsyslab::core::NpyArray::fromString(
        postprocess_detail::metadata(observed,method,label,details));
    return postprocess_detail::publish(root,observed,method,label,std::move(output));
}

}  // namespace dynsyslab::workflow
