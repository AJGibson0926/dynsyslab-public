#include "PostprocessInternal.hpp"
#include "ObservationSource.hpp"
#include "KromDatasetSupport.hpp"
#include "CoordinateInverse.hpp"
#include <dynsyslab/core/TimeSeriesOperations.hpp>
#include <dynsyslab/core/CoordinateTransform.hpp>
#include <dynsyslab/core/SindyDynamics.hpp>
#include <dynsyslab/core/DynamicsLibrary.hpp>
#include <dynsyslab/core/PathSimplification.hpp>
#include <utility>
#include <algorithm>
#include <iterator>
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
// q: predictor coordinates. numeric/reference: rates of the response
// coordinates (by default the predictors themselves) on the same clock.
struct Dataset {
    SourceData source;c::TimeSeries q,numeric;
    std::optional<c::TimeSeries> reference;
    bool cross{};
};
Dataset dataset(const SourceData& s) {
    if(datasetDynamics(s)!="unforced")
        throw std::invalid_argument("Autonomous dynamics SINDy requires confirmed unforced records");
    const auto series=c::TimeSeries::view(s.time,s.values.data(),s.samples,s.dimension);
    const auto rates=c::differentiateTimeSeries(series);
    Dataset d{s,c::cropSamples(series,2,s.samples-4),rates,{}};
    const auto a=c::readNpz(s.path);
    if(a.contains("coordinate_derivatives")) {
        const auto& v=a.at("coordinate_derivatives");
        if(v.shape!=std::vector<std::size_t>{s.samples,s.dimension})
            throw std::invalid_argument("Stored derivatives must align with coordinate samples");
        d.reference=c::cropSamples(c::TimeSeries(s.time,v.doubles(),s.dimension),2,s.samples-4);
    }
    return d;
}
// Rows of `series` at the given times (which must all be present).
c::TimeSeries rowsAt(const c::TimeSeries& series,const std::vector<double>& times) {
    std::vector<double> values;values.reserve(times.size()*series.componentCount());
    std::size_t k=0;
    for(const double t:times) {
        while(k<series.sampleCount()&&series.time(k)<t)++k;
        for(std::size_t j=0;j<series.componentCount();++j)values.push_back(series.value(k,j));
    }
    return {times,std::move(values),series.componentCount()};
}
// Predictors from one coordinate archive, response rates from another, on
// their common sample times (e.g. [e,e']' against (r,v) or xi).
Dataset paired(const SourceData& predictors,const SourceData& response) {
    const auto x=dataset(predictors),y=dataset(response);
    std::vector<double> common;
    std::set_intersection(x.q.timeData(),x.q.timeData()+x.q.sampleCount(),
        y.q.timeData(),y.q.timeData()+y.q.sampleCount(),std::back_inserter(common));
    if(common.size()<16)throw std::invalid_argument("Predictor and response coordinates share too few sample times");
    Dataset d{predictors,rowsAt(x.q,common),rowsAt(y.numeric,common),{},true};
    if(y.reference)d.reference=rowsAt(*y.reference,common);
    return d;
}
// Rows `kept` of several aligned series, archived under key+suffix, plus key+"time".
void archiveRows(c::NpzRecord& a,const std::string& key,const std::vector<std::pair<std::string,c::TimeSeries>>& series,
                 const std::vector<std::size_t>& kept) {
    std::vector<double> time;for(const auto k:kept)time.push_back(series.front().second.time(k));
    a.arrays[key+"time"]=c::NpyArray::fromDoubles({kept.size()},time);
    for(const auto& [suffix,s]:series) {
        std::vector<double> values;
        for(const auto k:kept)for(std::size_t j=0;j<s.componentCount();++j)values.push_back(s.value(k,j));
        a.arrays[key+suffix]=c::NpyArray::fromDoubles({kept.size(),s.componentCount()},values);
    }
}
// Error-bounded vertices of all the series together (core PathSimplification):
// collapses keep their detail, smooth stretches are thinned.
std::vector<std::size_t> rolloutVertices(const std::vector<std::pair<std::string,c::TimeSeries>>& series,std::size_t maximum) {
    const auto n=series.front().second.sampleCount();std::size_t width=0;
    for(const auto& item:series)width+=item.second.componentCount();
    std::vector<double> joined;joined.reserve(n*width);
    for(std::size_t k=0;k<n;++k)for(const auto& item:series)
        for(std::size_t j=0;j<item.second.componentCount();++j)joined.push_back(item.second.value(k,j));
    return c::simplifiedPathIndices(joined,n,width,maximum);
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
            c::dynamicsError(c::cropSamples(d.reference?*d.reference:d.numeric,first,last-first),
                c::cropSamples(prediction,first,last-first)).toJson()}});
    }
    result["regions"]=regions;
    json::Value::Array rollouts;std::size_t index=0;
    if(d.cross||!model.isVectorField()){result["rollouts"]="not applicable: the response coordinates differ from the predictors";return result;}
    for(const auto& start:settings.at("rollout_starts").array()) {
        const auto* begin=d.q.timeData();const auto* end=begin+d.q.sampleCount();
        const auto k=std::size_t(std::lower_bound(begin,end,start.number())-begin);
        // "full" runs to the end of the record; a longer duration is clipped to it.
        const auto& duration=settings.at("rollout_duration");
        const bool full=std::holds_alternative<std::string>(duration.storage())&&duration.string()=="full";
        if(k>=d.q.sampleCount())continue;
        auto count=full?d.q.sampleCount()-k:std::size_t(std::llround(duration.number(1)/d.q.deltaTime()))+1;
        count=std::min(count,d.q.sampleCount()-k);
        if(count<2)continue;
        json::Value::Object record{{"start",d.q.time(k)},{"duration",double(count-1)*d.q.deltaTime()}};
        try {
            // Integrate in segments (rollout_segment, default 0.5 tau), each continuing
            // from the previous end state, so a law that diverges keeps the valid
            // part of its trajectory instead of losing all of it.
            const std::size_t width=d.source.dimension;
            const auto full=c::cropSamples(d.q,k,count);
            const auto segment=std::max<std::size_t>(2,std::size_t(std::llround(
                settings.at("rollout_segment").number(0.5)/d.q.deltaTime())));
            std::vector<double> state(full.valuesData(),full.valuesData()+width),values=state;
            std::size_t valid=1,step=segment;std::string reason;
            while(valid<count) {
                const std::size_t a=valid-1,b=std::min(count-1,a+step);
                try {
                    c::IntegrationOptions options;options.initialTime=full.time(a);options.finalTime=full.time(b);
                    options.outputStep=d.q.deltaTime();
                    options.maximumInternalSteps=std::size_t(settings.at("rollout_max_steps").number(100000));
                    options.relativeTolerance=1e-8;options.absoluteTolerance=1e-10;
                    const auto solution=c::integrateSindyDynamics(model,state,options);
                    if(solution.sampleCount<b-a+1)throw std::runtime_error("Rollout output grid differs from requested horizon");
                    // Compare on the observed timestamps (solver roundoff need not be bitwise identical).
                    for(std::size_t row=1;row<=b-a;++row) {
                        if(std::abs(solution.time[row]-full.time(a+row))>1e-10*std::max(1.,std::abs(full.time(a+row))))
                            throw std::runtime_error("Rollout timestamps do not match validation samples");
                        values.insert(values.end(),solution.values.begin()+std::ptrdiff_t(row*width),
                                      solution.values.begin()+std::ptrdiff_t((row+1)*width));
                    }
                    state.assign(solution.values.begin()+std::ptrdiff_t((b-a)*width),solution.values.begin()+std::ptrdiff_t((b-a+1)*width));
                    valid=b+1;step=std::min(segment,2*step);
                } catch(const std::exception& e) {
                    // Halve the segment to keep everything up to the blow-up itself.
                    reason=e.what();
                    if(step<=2)break;
                    step=std::max<std::size_t>(2,step/2);
                }
            }
            if(valid<2)throw std::runtime_error(reason.empty()?"rollout produced no samples":reason);
            const auto actual=c::cropSamples(d.q,k,valid);
            const auto predicted=c::TimeSeries({actual.timeData(),actual.timeData()+valid},std::move(values),width);
            record["error"]=c::dynamicsError(actual,predicted).toJson();
            record["status"]=valid==count?"complete":"diverged";
            if(valid<count){record["valid_until"]=actual.time(valid-1);record["reason"]=reason;}
            const auto name=prefix+"rollout_"+std::to_string(index)+"_";
            std::vector<std::pair<std::string,c::TimeSeries>> stored{{"actual",actual},{"predicted",predicted}};
            // rollout_map: show a state-space rollout in other coordinates too, e.g.
            // echo_jet2 maps a free (r,v) rollout and its reference to (e,e').
            if(const auto mapName=settings.at("rollout_map").string();!mapName.empty()) {
                const auto catalog=c::createSystem(systemConfiguration(d.source))->observableCatalog();
                const auto map=catalog.coordinatePresets().count(mapName)?catalog.coordinatePreset(mapName):catalog.at(mapName);
                stored.emplace_back("mapped_actual",c::transformCoordinates(actual,map).coordinates);
                stored.emplace_back("mapped_predicted",c::transformCoordinates(predicted,map).coordinates);
                record["mapped_error"]=c::dynamicsError(stored[2].second,stored[3].second).toJson();
                record["map"]=mapName;
                const auto& names=map.descriptor().outputNames;
                record["map_labels"]=names.size()==2&&mapName=="echo_jet2"?json::Value::Array{"e","e'"}
                    :json::parse(jsonStringArray(names));
            } else if(const auto inversePath=settings.at("rollout_inverse").string();!inversePath.empty()) {
                // The coordinates' retained inverse, e.g. xi back to (e,e'): both
                // the observed and the free rollout go through the same map.
                const auto inverse=loadCoordinateInverse(std::filesystem::u8path(inversePath));
                stored.emplace_back("mapped_actual",inverse.apply(actual));
                stored.emplace_back("mapped_predicted",inverse.apply(predicted));
                record["mapped_error"]=c::dynamicsError(stored[2].second,stored[3].second).toJson();
                record["map"]="inverse ("+inverse.kind()+")";
                record["map_labels"]=json::parse(jsonStringArray(inverse.targetLabels()));
            }
            // Stored for display with error-bounded vertices (rollout_vertices),
            // not every Nth sample; errors above use every sample.
            archiveRows(a,name,stored,rolloutVertices(stored,std::size_t(settings.at("rollout_vertices").number(60000))));
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
    // response_source: fit the rates of another coordinate archive's coordinates
    // (w' = g(z)) instead of the predictors' own rates (z' = f(z)).
    const auto responsePath=settings.at("response_source").string();
    const auto data=responsePath.empty()?dataset(source):paired(source,[&]{
        auto r=loadSource(std::filesystem::u8path(responsePath));requireSameObservationSystem(source,r);return r;}());
    const auto dimension=data.numeric.componentCount();
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
    archive.arrays["coefficients"]=c::NpyArray::fromDoubles({library.outputDimension(),model.outputs()},model.coefficients());
    auto meta=json::parse(metadata(source,"sindy-dynamics",label,""));
    const auto sourceMeta=json::parse(source.metadataJson);
    for(const char* key:{"coordinate_names","coordinate_labels","coordinate_transform_kind","coordinate_transform_identity",
        "coordinate_transform_archive","coordinate_definition","derivative_evidence"})
        if(!sourceMeta.at(key).isNull())meta.object()[key]=sourceMeta.at(key);
    meta.object()["settings"]=settings;meta.object()["coordinate_dimension"]=double(source.dimension);
    meta.object()["response_dimension"]=double(dimension);
    if(!responsePath.empty()) {
        const auto r=json::parse(c::readNpz(std::filesystem::u8path(responsePath)).metadataJson());
        meta.object()["response_source"]=responsePath;
        meta.object()["response_coordinate_names"]=r.at("coordinate_names");
        meta.object()["response_coordinate_labels"]=r.at("coordinate_labels");
    }
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
    const auto& responses=settings.at("validation_response_sources").array();
    if(!responsePath.empty()&&responses.size()!=settings.at("validation_sources").array().size())
        throw std::invalid_argument("validation_response_sources must pair with validation_sources");
    for(const auto& path:settings.at("validation_sources").array()) {
        const auto s=loadSource(std::filesystem::u8path(path.string()));
        requireSameObservationSystem(source,s);
        const auto first=json::parse(source.metadataJson),next=json::parse(s.metadataJson);
        for(const char* key:{"coordinate_transform_kind","coordinate_transform_identity","coordinate_definition"})
            if(json::dump(first.at(key))!=json::dump(next.at(key)))
                throw std::invalid_argument("Validation coordinates must use the frozen training transform");
        if(s.dimension!=source.dimension)throw std::invalid_argument("Validation coordinate dimension mismatch");
        const auto validation=responsePath.empty()?dataset(s)
            :paired(s,loadSource(std::filesystem::u8path(responses[index].string())));
        validations.push_back(evaluate(model,validation,archive,"validation_"+std::to_string(index++)+"_",settings));
    }
    meta.object()["validation"]=validations;
    meta.object()["model_type"]=model.isVectorField()&&responsePath.empty()?"autonomous-vector-field":"cross-coordinate-rate-law";
    meta.object()["derivative_method"]=derivative;
    archive.arrays["metadata_json"]=c::NpyArray::fromString(json::dump(meta));
    return publish(root,source,"sindy-dynamics",label,std::move(archive));
}
}
