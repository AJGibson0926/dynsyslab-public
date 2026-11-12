#include "PostprocessInternal.hpp"
#include "ObservationSource.hpp"
#include "KromDatasetSupport.hpp"
#include <dynsyslab/core/CoordinateTransform.hpp>
#include <dynsyslab/core/SindyDynamics.hpp>
#include <dynsyslab/core/DynamicsLibrary.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace {
namespace c=dynsyslab::core;
std::size_t positive(const json::Value& settings,const char* key,std::size_t fallback) {
    const auto v=settings.at(key).number(double(fallback));
    if(!std::isfinite(v)||v<1||v!=std::floor(v)||v>1e8)
        throw std::invalid_argument(std::string("Invalid positive integer: ")+key);
    return std::size_t(v);
}
c::TimeSeries sampled(const c::TimeSeries& s,std::size_t stride) {
    std::vector<double> time,values;
    for(std::size_t k=0;k<s.sampleCount();k+=stride) {
        time.push_back(s.time(k));
        for(std::size_t j=0;j<s.componentCount();++j)values.push_back(s.value(k,j));
    }
    return {time,values,s.componentCount()};
}
struct Dataset {
    SourceData source;c::TimeSeries q,numeric;
    std::optional<c::TimeSeries> reference;
};
Dataset dataset(const SourceData& s) {
    if(datasetDynamics(s)!="unforced")
        throw std::invalid_argument("Autonomous dynamics SINDy requires confirmed unforced records");
    const auto series=c::TimeSeries::view(s.time,s.values.data(),s.samples,s.dimension);
    const auto rates=c::differentiateTimeSeries(series);
    Dataset d{s,c::trimTimeSeries(series,2,s.samples-4),rates,{}};
    const auto a=c::readNpz(s.path);
    if(a.contains("coordinate_derivatives")) {
        const auto& v=a.at("coordinate_derivatives");
        if(v.shape!=std::vector<std::size_t>{s.samples,s.dimension})
            throw std::invalid_argument("Stored derivatives must align with coordinate samples");
        d.reference=c::trimTimeSeries(c::TimeSeries(s.time,v.doubles(),s.dimension),2,s.samples-4);
    }
    return d;
}
void archiveSeries(c::NpzRecord& a,const std::string& key,const c::TimeSeries& s,std::size_t stride) {
    auto thin=sampled(s,stride);
    a.arrays[key]=c::NpyArray::fromDoubles({thin.sampleCount(),thin.componentCount()},
        {thin.valuesData(),thin.valuesData()+thin.sampleCount()*thin.componentCount()});
}
json::Value evaluate(const c::SindyDynamicsModel& model,const Dataset& d,
    c::NpzRecord& a,const std::string& prefix,const json::Value& settings) {
    const auto prediction=model.evaluate(d.q);const auto stride=positive(settings,"diagnostic_stride",1);
    const auto thin=sampled(d.q,stride);
    a.arrays[prefix+"time"]=c::NpyArray::fromDoubles({thin.sampleCount()},
        {thin.timeData(),thin.timeData()+thin.sampleCount()});
    archiveSeries(a,prefix+"coordinates",d.q,stride);
    archiveSeries(a,prefix+"predicted_derivatives",prediction,stride);
    archiveSeries(a,prefix+"numerical_derivatives",d.numeric,stride);
    json::Value::Object result{{"source_archive",std::filesystem::absolute(d.source.path).generic_u8string()},
        {"against_numerical",c::dynamicsError(d.numeric,prediction).toJson()}};
    if(d.reference) {
        archiveSeries(a,prefix+"reference_derivatives",*d.reference,stride);
        result["against_reference"]=c::dynamicsError(*d.reference,prediction).toJson();
        result["differentiation_error"]=c::dynamicsError(*d.reference,d.numeric).toJson();
    }
    json::Value::Array regions;
    for(const auto& interval:settings.at("evaluation_intervals").array()) {
        if(interval.array().size()!=2)throw std::invalid_argument("Evaluation interval needs start and stop");
        const double start=interval.array()[0].number(),stop=interval.array()[1].number();
        const auto* begin=d.q.timeData();const auto* end=begin+d.q.sampleCount();
        const auto first=std::lower_bound(begin,end,start)-begin,last=std::upper_bound(begin,end,stop)-begin;
        if(last<=first)continue;
        regions.emplace_back(json::Value::Object{{"start",start},{"stop",stop},{"error",
            c::dynamicsError(c::trimTimeSeries(d.reference?*d.reference:d.numeric,first,last-first),
                c::trimTimeSeries(prediction,first,last-first)).toJson()}});
    }
    result["regions"]=regions;
    json::Value::Array rollouts;std::size_t index=0;
    for(const auto& start:settings.at("rollout_starts").array()) {
        const auto* begin=d.q.timeData();const auto* end=begin+d.q.sampleCount();
        const auto k=std::size_t(std::lower_bound(begin,end,start.number())-begin);
        const auto count=std::size_t(std::llround(settings.at("rollout_duration").number(1)/d.q.deltaTime()))+1;
        if(k>=d.q.sampleCount()||count>d.q.sampleCount()-k||count<2)continue;
        json::Value::Object record{{"start",d.q.time(k)},{"duration",double(count-1)*d.q.deltaTime()}};
        try {
            c::IntegrationOptions options;options.initialTime=d.q.time(k);options.finalTime=d.q.time(k+count-1);
            options.outputStep=d.q.deltaTime();options.maximumInternalSteps=100000;
            options.relativeTolerance=1e-8;options.absoluteTolerance=1e-10;
            std::vector<double> initial(d.q.valuesData()+k*d.source.dimension,d.q.valuesData()+(k+1)*d.source.dimension);
            const auto solution=c::integrateSindyDynamics(model,initial,options);
            // Compare using the observed timestamps (solver roundoff need not be bitwise identical).
            const auto actual=c::trimTimeSeries(d.q,k,count);
            if(solution.sampleCount<count||solution.sampleCount>count+1)
                throw std::runtime_error("Rollout output grid differs from requested horizon");
            for(std::size_t row=0;row<count;++row)
                if(std::abs(solution.time[row]-actual.time(row))>1e-10*std::max(1.,std::abs(actual.time(row))))
                    throw std::runtime_error("Rollout timestamps do not match validation samples");
            // The standard solver's ceil grid may include one duplicate endpoint
            // from floating-point rounding. Preserve one copy of each observed tick.
            const auto predicted=c::TimeSeries({actual.timeData(),actual.timeData()+count},
                {solution.values.data(),solution.values.data()+count*d.source.dimension},d.source.dimension);
            record["error"]=c::dynamicsError(actual,predicted).toJson();record["status"]="complete";
            const auto name=prefix+"rollout_"+std::to_string(index)+"_";
            archiveSeries(a,name+"actual",actual,stride);archiveSeries(a,name+"predicted",predicted,stride);
            const auto ts=sampled(actual,stride);
            a.arrays[name+"time"]=c::NpyArray::fromDoubles({ts.sampleCount()},{ts.timeData(),ts.timeData()+ts.sampleCount()});
        } catch(const std::exception& e){record["status"]="failed";record["reason"]=e.what();}
        rollouts.emplace_back(record);++index;
    }
    result["rollouts"]=rollouts;return result;
}
}
std::filesystem::path runDynamicsValidation(const SourceData& source,const json::Value& settings,
    const std::filesystem::path& root,const std::string& label) {
    const auto path=std::filesystem::u8path(settings.at("model_archive").string());
    const auto saved=c::readNpz(path);
    const auto model=c::SindyDynamicsModel::fromDefinition(json::parse(saved.at("dynamics_model_json").string()));
    const SourceData provenance{path,{},{},0,model.dimension(),{},0,saved.metadataJson()};
    requireSameObservationSystem(source,provenance);
    const auto sourceMeta=json::parse(source.metadataJson),modelMeta=json::parse(saved.metadataJson());
    for(const char* key:{"coordinate_transform_kind","coordinate_transform_identity","coordinate_definition"})
        if(json::dump(sourceMeta.at(key))!=json::dump(modelMeta.at(key)))
            throw std::invalid_argument("Model validation requires the original coordinate transform");
    const auto data=dataset(source);
    c::NpzRecord archive;
    auto meta=json::parse(metadata(source,"dynamics-validation",label,""));
    meta.object()["model_archive"]=std::filesystem::absolute(path).generic_u8string();
    meta.object()["coordinate_dimension"]=double(source.dimension);
    meta.object()["validation"]=evaluate(model,data,archive,"validation_",settings);
    if(!settings.at("compare_model_archive").isNull()) {
        const auto comparison=c::readNpz(std::filesystem::u8path(settings.at("compare_model_archive").string()));
        const auto other=c::SindyDynamicsModel::fromDefinition(json::parse(comparison.at("dynamics_model_json").string()));
        const auto otherMeta=json::parse(comparison.metadataJson());
        const SourceData otherProvenance{{},{},{},0,other.dimension(),{},0,comparison.metadataJson()};
        requireSameObservationSystem(source,otherProvenance);
        for(const char* key:{"coordinate_transform_kind","coordinate_transform_identity","coordinate_definition"})
            if(json::dump(sourceMeta.at(key))!=json::dump(otherMeta.at(key)))
                throw std::invalid_argument("Comparison model uses different coordinates");
        meta.object()["model_comparison"]=c::dynamicsError(other.evaluate(data.q),model.evaluate(data.q)).toJson();
        meta.object()["comparison_model_archive"]=settings.at("compare_model_archive");
    }
    archive.arrays["metadata_json"]=c::NpyArray::fromString(json::dump(meta));
    return publish(root,source,"dynamics-validation",label,std::move(archive));
}
std::filesystem::path runSindyDynamics(const SourceData& source,const json::Value& settings,
    const std::filesystem::path& root,const std::string& label) {
    const auto data=dataset(source);const auto dimension=source.dimension;
    const auto derivative=settings.at("derivatives").string("finite_difference");
    if(derivative!="finite_difference"&&derivative!="stored")
        throw std::invalid_argument("SINDy dynamics derivatives must be finite_difference or stored");
    if(derivative=="stored"&&!data.reference)throw std::invalid_argument("Missing stored derivative reference");
    const auto& rates=derivative=="stored"?*data.reference:data.numeric;
    const auto stride=positive(settings,"fit_stride",1);
    auto q=sampled(data.q,stride),dq=sampled(rates,stride);
    const auto family=settings.at("library_family").string("broad");
    std::vector<std::string> exclusions;
    const auto library=[&] {
        if(!settings.at("library").isNull())return c::Observable::fromExpressionDefinition(settings.at("library"));
        if(family!="broad"&&family!="polynomial")throw std::invalid_argument("library_family must be broad or polynomial (or supply library expressions)");
        c::DynamicsLibraryOptions o;o.polynomialDegree=positive(settings,"polynomial_degree",4);
        o.normalize=settings.at("normalize").boolean(true);const auto& f=settings.at("library_families");
        o.trigonometric=f.at("trigonometric").boolean(true);o.reciprocal=f.at("reciprocal").boolean(true);
        o.logarithmic=f.at("logarithmic").boolean(true);o.exponential=f.at("exponential").boolean(true);
        o.hyperbolic=f.at("hyperbolic").boolean(true);o.mixedTerms=f.at("mixed_terms").boolean(true);
        if(family=="polynomial")o.trigonometric=o.reciprocal=o.logarithmic=o.exponential=o.hyperbolic=o.mixedTerms=false;
        o.inverseDegree=positive(settings,"inverse_degree",5);
        auto built=c::buildDynamicsLibrary(data.q,o);exclusions=built.excluded;return built.observable;
    }();
    c::SindyDynamicsOptions options;options.threshold=settings.at("threshold").number(.001);
    options.maximumIterations=positive(settings,"max_iterations",20);
    options.singularValueTolerance=settings.at("svd_tolerance").number(1e-10);
    c::SindyDynamicsFitDiagnostics fitDiagnostics;
    const auto model=c::fitSindyDynamics(q,dq,library,options,&fitDiagnostics);
    c::NpzRecord archive;
    archive.arrays["dynamics_model_json"]=c::NpyArray::fromString(json::dump(model.definition()));
    archive.arrays["coefficients"]=c::NpyArray::fromDoubles({library.outputDimension(),dimension},model.coefficients());
    auto meta=json::parse(metadata(source,"sindy-dynamics",label,""));
    const auto sourceMeta=json::parse(source.metadataJson);
    for(const char* key:{"coordinate_names","coordinate_transform_kind","coordinate_transform_identity",
        "coordinate_transform_archive","coordinate_definition","derivative_evidence"})
        if(!sourceMeta.at(key).isNull())meta.object()[key]=sourceMeta.at(key);
    meta.object()["settings"]=settings;meta.object()["coordinate_dimension"]=double(dimension);
    meta.object()["fit_samples"]=double(q.sampleCount());
    meta.object()["candidate_count"]=double(library.outputDimension());
    meta.object()["library_family"]=settings.at("library").isNull()?family:"custom";
    meta.object()["excluded_candidates"]=json::parse(jsonStringArray(exclusions));
    json::Value::Array support;
    for(std::size_t j=0;j<dimension;++j) {
        json::Value::Array active;
        for(std::size_t k=0;k<library.outputDimension();++k)if(model.coefficients()[k*dimension+j]!=0)
            active.emplace_back(library.descriptor().outputNames.size()>k?
                library.descriptor().outputNames[k]:"feature "+std::to_string(k));
        support.emplace_back(active);
    }
    meta.object()["selected_terms"]=support;
    meta.object()["sparse_regression_steps"]=fitDiagnostics.toJson();
    meta.object()["training"]=evaluate(model,data,archive,"training_",settings);
    json::Value::Array validations;std::size_t index=0;
    for(const auto& path:settings.at("validation_sources").array()) {
        const auto s=loadSource(std::filesystem::u8path(path.string()));
        requireSameObservationSystem(source,s);
        const auto first=json::parse(source.metadataJson),next=json::parse(s.metadataJson);
        for(const char* key:{"coordinate_transform_kind","coordinate_transform_identity","coordinate_definition"})
            if(json::dump(first.at(key))!=json::dump(next.at(key)))
                throw std::invalid_argument("Validation coordinates must use the frozen training transform");
        if(s.dimension!=dimension)throw std::invalid_argument("Validation coordinate dimension mismatch");
        validations.push_back(evaluate(model,dataset(s),archive,"validation_"+std::to_string(index++)+"_",settings));
    }
    meta.object()["validation"]=validations;
    meta.object()["model_type"]="autonomous-vector-field";
    meta.object()["derivative_method"]=derivative;
    archive.arrays["metadata_json"]=c::NpyArray::fromString(json::dump(meta));
    return publish(root,source,"sindy-dynamics",label,std::move(archive));
}
}
