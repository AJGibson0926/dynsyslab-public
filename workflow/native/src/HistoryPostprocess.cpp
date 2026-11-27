#include "MeasurementMap.hpp"
#include <dynsyslab/workflow/HistoryPostprocess.hpp>
#include "ObservationSource.hpp"
#include <dynsyslab/core/CoordinateTransform.hpp>
#include "KromDatasetSupport.hpp"
#include <dynsyslab/core/HistoryKoopman.hpp>
#include <dynsyslab/core/HistoryObservationArchive.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace dynsyslab::workflow {
namespace {
namespace c=dynsyslab::core;
namespace p=postprocess_detail;
std::size_t count(const json::Value& s,const char* key,std::size_t fallback) {
    const double n=s.at(key).number(double(fallback));
    if(!std::isfinite(n)||n<1||n!=std::floor(n)||n>double(std::numeric_limits<std::size_t>::max()/2))
        throw std::invalid_argument(std::string("Invalid positive count: ")+key);
    return static_cast<std::size_t>(n);
}
p::SourceData sourceData(const std::filesystem::path& path,const json::Value& settings) {
    auto s=p::loadObservationSource(path,settings.at("observation_array").string("echo"),settings.at("observable").string());
    if(s.dimension()!=1)throw std::invalid_argument("History transform requires scalar observations");
    return s;
}
std::vector<double> decimate(const std::vector<double>& x,std::size_t stride) {
    std::vector<double> out;for(std::size_t k=0;k<x.size();k+=stride)out.push_back(x[k]);return out;
}
std::string measurementName(const json::Value& settings) {
    return settings.at("observable").string(settings.at("observation_array").string("echo"));
}
std::string observationUnit(const p::SourceData& source,const json::Value& settings,
    const char* setting,const char* metadataKey,const char* fallback) {
    const auto meta=source.metadataJson.empty()?json::Value{}:json::parse(source.metadataJson);
    const auto declared=meta.at(metadataKey).string();
    const auto selected=settings.at(setting).string(declared.empty()?fallback:declared);
    if(!declared.empty()&&declared!=selected)
        throw std::invalid_argument(std::string("Observation unit override conflicts with source: ")+setting);
    return selected;
}
}
std::filesystem::path runHistoryBasis(const std::vector<std::filesystem::path>& paths,
    const json::Value& settings,const std::filesystem::path& root,const std::string& label) {
    if(paths.empty())throw std::invalid_argument("History fitting needs source observations");
    const auto first=sourceData(paths.front(),settings);
    const auto stride=count(settings,"decimation",1);
    const double dt=p::uniformDeltaTime(first)*double(stride);
    std::vector<std::vector<double>> observations;
    std::vector<c::KoopmanDatasetReference> references;
    for(const auto& path:paths) {
        const auto source=sourceData(path,settings);
        if(p::datasetDynamics(source)!="unforced")throw std::invalid_argument("History basis training requires confirmed unforced observations");
        p::requireSameObservationSystem(first,source);
        if(std::abs(p::uniformDeltaTime(source)*double(stride)-dt)>1e-8*dt)
            throw std::invalid_argument("Pooled history clocks differ");
        observations.push_back(decimate(source.values(),stride));
        references.push_back({std::filesystem::absolute(path).generic_u8string(),"training","unforced",{}});
    }
    c::PooledHankelOptions options;options.delays=count(settings,"delays",100);
    options.rank=count(settings,"hankel_rank",4);options.oversampling=count(settings,"hankel_oversampling",4);
    options.subspaceIterations=count(settings,"hankel_subspace_iterations",8);options.seed=count(settings,"seed",20261002);
    auto basis=c::fitPooledHankelBasis(observations,options);
    const auto rank=count(settings,"coordinate_rank",options.rank);
    if(rank>basis.rank)throw std::invalid_argument("Coordinate rank exceeds the fitted basis");
    std::vector<std::size_t> selection;for(std::size_t j=0;j<rank;++j)selection.push_back(j);
    const c::FrozenHistoryTransform transform(std::move(basis),dt,selection,
        measurementName(settings),
        observationUnit(first,settings,"time_units","measurement_time_units","model-time"),
        observationUnit(first,settings,"observation_units","measurement_units","unspecified"));
    c::NpzRecord a;c::writeFrozenHistoryTransform(a,transform);
    auto meta=json::parse(p::metadata(first,"history-basis",label,""));
    meta.object()["dataset_references"]=c::koopmanDatasetReferencesJson(references);
    meta.object()["history_decimation"]=double(stride);
    meta.object()["transform_identity"]=transform.identity();
    a.arrays["metadata_json"]=c::NpyArray::fromString(json::dump(meta));
    return p::publishObject(root, first, dynsyslab::core::object_kind::coordinateMap, "history-basis",std::move(a));
}
std::filesystem::path runHistoryProjection(const std::filesystem::path& path,
    const json::Value& settings,const std::filesystem::path& root,const std::string& label) {
    const auto transformPath=std::filesystem::u8path(settings.at("transform").string());
    if(transformPath.empty())throw std::invalid_argument("History projection requires a frozen transform archive");
    const auto transform=postprocess_detail::readMeasurementMap(transformPath);
    const auto source=sourceData(path,settings);
    if(measurementName(settings)!=transform.measurement())
        throw std::invalid_argument("Projection measurement differs from the frozen transform");
    if(observationUnit(source,settings,"time_units","measurement_time_units",transform.timeUnits().c_str())!=transform.timeUnits()
        ||observationUnit(source,settings,"observation_units","measurement_units",transform.observationUnits().c_str())!=transform.observationUnits())
        throw std::invalid_argument("Projection units differ from the frozen transform; convert observations explicitly");
    const auto stride=count(settings,"decimation",1);
    if(std::abs(p::uniformDeltaTime(source)*double(stride)-transform.sampleInterval())>1e-8*transform.sampleInterval())
        throw std::invalid_argument("Projection observation clock differs from the frozen transform");
    const auto observations=decimate(source.values(),stride),time=decimate(source.time(),stride);
    const auto projected=c::transformHistoryCoordinates(c::TimeSeries(time,observations,1),transform);
    const auto& selection=transform.retainedCoordinates();
    const auto& q=projected.coordinates;
    const auto windows=q.sampleCount();
    c::NpzRecord a;
    a.arrays["time"]=c::NpyArray::fromDoubles({windows},{q.timeData(),q.timeData()+windows});
    a.arrays["values"]=c::NpyArray::fromDoubles({windows,selection.size()},
        {q.valuesData(),q.valuesData()+windows*selection.size()});
    {   // The physical states and the measurement at each window's newest sample, for audits that
        // compare coordinates with the states they stand for (orbit_geometry).
        const auto raw=c::readNpz(path);
        std::string stateKey=raw.contains("source_states")?"source_states":"";
        if(stateKey.empty()&&raw.contains("values")&&settings.at("observation_array").string("echo")!="values")stateKey="values";
        std::size_t offset=0;
        while(offset<time.size()&&windows&&time[offset]<q.timeData()[0]-1e-9*std::max(1.0,std::abs(q.timeData()[0])))++offset;
        std::vector<double> measured;
        for(std::size_t k=0;k<windows&&offset+k<observations.size();++k)measured.push_back(observations[offset+k]);
        if(measured.size()==windows)a.arrays["measurement_values"]=c::NpyArray::fromDoubles({windows},measured);
        if(!stateKey.empty()) {
            const auto& states=raw.at(stateKey);
            const auto width=states.shape.size()==2?states.shape[1]:0;
            if(width&&states.shape[0]==source.samples()&&(offset+windows-1)*stride<states.shape[0]) {
                const auto all=states.doubles();
                std::vector<double> carried;
                for(std::size_t k=0;k<windows;++k){const auto row=(offset+k)*stride;carried.insert(carried.end(),all.begin()+std::ptrdiff_t(row*width),all.begin()+std::ptrdiff_t((row+1)*width));}
                a.arrays["source_states"]=c::NpyArray::fromDoubles({windows,width},carried);
            }
        }
    }
    auto meta=json::parse(p::metadata(source,"history-projection",label,""));
    meta.object()["coordinate_provider"]="causal-history";
    meta.object()["coordinate_transform_identity"]=transform.identity();
    meta.object()["coordinate_transform_archive"]=std::filesystem::absolute(transformPath).generic_u8string();
    meta.object()["value_space"]="history-coordinates";
    const auto dynamics=p::datasetDynamics(source);
    meta.object()["forcing_type"]=dynamics=="unforced"?"free":dynamics;
    meta.object()["controller_kind"]=dynamics=="unforced"?"none":"unknown";
    meta.object()["warmup_samples"]=double(transform.basis().rows-1);
    meta.object()["delta_t"]=transform.sampleInterval();
    // Original controls remain in the source archive. No implicit input averaging.
    meta.object()["input_history_source"]=std::filesystem::absolute(path).generic_u8string();
    a.arrays["metadata_json"]=c::NpyArray::fromString(json::dump(meta));
    return p::publishObject(root, source, dynsyslab::core::object_kind::timeSeries, "history-projection",std::move(a));
}
}
