#include "PostprocessInternal.hpp"
#include "ObservationSource.hpp"
#include "KromDatasetSupport.hpp"
#include <dynsyslab/core/AcousticReceiver.hpp>
#include <dynsyslab/core/CoordinateTransform.hpp>
#include <dynsyslab/core/HistoryObservationArchive.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace {
core::Observable receiverMap(const SourceData& source, const json::Value& settings, json::Value& meta) {
    const auto& spec=settings.at("receiver");
    const auto& pos=spec.at("position_m").array();
    if(pos.size()!=3)throw std::invalid_argument("receiver.position_m needs three offsets in metres");
    const double missing=std::numeric_limits<double>::quiet_NaN();
    const core::AcousticReceiver receiver{{pos[0].number(missing),pos[1].number(missing),pos[2].number(missing)},
        spec.at("sound_speed_m_s").number(missing)};
    const auto system=core::createSystem(systemConfiguration(source));
    const auto catalog=core::rpeReceiverPressures(*system,receiver);
    const auto name=settings.at("observable").string();
    const auto& map=catalog.at(name);
    const double range=receiver.rangeMetres();
    const double radius0=system->configuration().parameters.at("equilibrium_radius_m");
    if(source.dimension!=2)throw std::invalid_argument("Receiver pressures require RPE radius and velocity");
    for(std::size_t row=0;row<source.samples;++row) {
        const double radius=radius0*source.values[row*2];
        if(!(radius>0 && radius<range))
            throw std::invalid_argument("Receiver must be outside the bubble at every recorded state");
    }
    auto saved=spec;
    saved.object()["term"]=name;
    saved.object()["range_m"]=range;
    saved.object()["travel_time_s"]=receiver.travelTimeSeconds();
    saved.object()["clock"]="emission-time";
    meta.object()["acoustic_receiver"]=std::move(saved);
    return map;
}
}
std::filesystem::path runCoordinateTransform(const SourceData& source,const json::Value& settings,
    const std::filesystem::path& root,const std::string& label) {
    namespace c=dynsyslab::core;
    const auto kind=settings.at("kind").string("identity");
    const bool modelRates=settings.at("system_derivatives").boolean(false);
    const bool integral=kind=="window_integral";
    const bool mapped=!settings.at("observable").isNull()||!settings.at("definition").isNull();
    if(integral&&modelRates)
        throw std::invalid_argument("Window integration produces coordinate values only");
    if(integral&&mapped&&!settings.at("measurement_array").isNull())
        throw std::invalid_argument("Choose either a measured array or an analytic observable for window integration");
    const auto original=source.metadataJson.empty()?json::Value{}:json::parse(source.metadataJson);
    const bool receiver=!settings.at("receiver").isNull();
    if(receiver&&((kind!="observable"&&!integral)||settings.at("observable").string().empty()||
       !settings.at("definition").isNull()||!settings.at("measurement_array").isNull()||
       !original.at("coordinate_provider").isNull()))
        throw std::invalid_argument("A receiver requires one named acoustic term and a physical state source");
    if(integral&&mapped&&!original.at("coordinate_provider").isNull())
        throw std::invalid_argument("Use the scalar measurement array on a transformed trajectory");
    if(modelRates&&datasetDynamics(source)!="unforced")
        throw std::invalid_argument("Autonomous derivative reference requires confirmed unforced data");
    if(modelRates&&!json::parse(source.metadataJson).at("coordinate_provider").isNull())
        throw std::invalid_argument("System derivatives require physical-state source data");
    const auto q=c::TimeSeries::view(source.time,source.values.data(),source.samples,source.dimension);
    std::optional<c::TimeSeries> rates;
    if(modelRates)rates=c::TimeSeries(source.time,exactVelocities(source,false),source.dimension);
    auto result=c::CoordinateSeries{q,rates};
    auto meta=json::parse(metadata(source,"coordinate-transform",label,""));
    meta.object()["coordinate_provider"]="measured";
    meta.object()["value_space"]="observation-coordinates";
    meta.object()["forcing_type"]=datasetDynamics(source)=="unforced"?"free":"unknown";
    meta.object()["derivative_evidence"]=modelRates?"system-rhs/chain-rule":"none";
    if(!settings.at("names").isNull())meta.object()["coordinate_names"]=settings.at("names");
    if(kind=="observable"||(integral&&mapped)) {
        auto map=[&] {
            if(receiver)return receiverMap(source,settings,meta);
            if(!settings.at("definition").isNull())return c::Observable::fromExpressionDefinition(settings.at("definition"));
            const auto catalog=c::createSystem(systemConfiguration(source))->observableCatalog();
            const auto name=settings.at("observable").string();
            return catalog.coordinatePresets().count(name)?catalog.coordinatePreset(name):catalog.at(name);
        }();
        result=c::transformCoordinates(q,map,rates?&*rates:nullptr);
        meta.object()["coordinate_definition"]=map.expressionDefinition();
    }
    if(integral) {
        const double duration=settings.at("window_duration").number();
        result.coordinates=c::windowIntegralCoordinates(result.coordinates,duration);
        json::Value::Object definition{{"kind","window_integral"},{"window_duration",duration},
            {"quadrature","piecewise-linear"},{"alignment","trailing-complete"}};
        if(mapped) definition["observable"]=meta.at("coordinate_definition");
        else {
            definition["measurement_array"]=settings.at("measurement_array").string("values");
            for(const auto key:{"coordinate_definition","coordinate_transform_identity","coordinate_names"})
                if(!original.at(key).isNull())definition[key]=original.at(key);
        }
        for(const auto key:{"measurement_units","measurement_time_units","model_time_units"})
            if(!original.at(key).isNull())definition[key]=original.at(key);
        meta.object()["coordinate_provider"]="causal-history";
        meta.object()["value_space"]="history-coordinates";
        meta.object()["coordinate_definition"]=definition;
        meta.object()["coordinate_transform_identity"]=json::dump(definition);
        meta.object()["history_span"]=duration;
        if(settings.at("names").isNull())
            meta.object()["coordinate_names"]=json::Value::Array{"signal","window_integral"};
    } else if(kind=="history") {
        const auto path=std::filesystem::u8path(settings.at("transform_archive").string());
        const auto transform=c::loadFrozenHistoryTransform(path);
        const auto observable=settings.at("observable").string(transform.measurement());
        if(observable!=transform.measurement()) {
            const auto system=c::createSystem(systemConfiguration(source));
            if(system->observableCatalog().at(observable).descriptor().name!=
                system->observableCatalog().at(transform.measurement()).descriptor().name)
                throw std::invalid_argument("Measurement differs from the frozen history transform");
        }
        for(const auto& pair:std::vector<std::pair<std::string,std::string>>{
            {"measurement_time_units",transform.timeUnits()},{"measurement_units",transform.observationUnits()}})
            if(!original.at(pair.first).isNull()&&original.at(pair.first).string()!=pair.second)
                throw std::invalid_argument("History measurement units differ from the frozen transform");
        const auto measured=loadObservationSource(source.path,settings.at("measurement_array").string("echo_values"),observable);
        auto measurement=c::TimeSeries::view(measured.time,measured.values.data(),measured.samples,measured.dimension);
        std::optional<c::TimeSeries> measurementRates;
        if(modelRates) {
            const auto map=c::createSystem(systemConfiguration(source))->observableCatalog().at(observable);
            // Exact reference is meaningful only when samples are this analytic observable.
            const auto analytic=c::transformCoordinates(q,map,&*rates);
            if(analytic.coordinates.componentCount()!=1||measurement.componentCount()!=1)
                throw std::invalid_argument("History transform requires one scalar observable");
            for(std::size_t k=0;k<measurement.sampleCount();++k)
                if(std::abs(measurement.value(k,0)-analytic.coordinates.value(k,0))>
                    1e-10*std::max(1.,std::abs(analytic.coordinates.value(k,0))))
                    throw std::invalid_argument("Measured signal differs from requested analytic derivative reference");
            measurementRates=analytic.derivatives;
        }
        result=c::transformHistoryCoordinates(measurement,transform,measurementRates?&*measurementRates:nullptr);
        meta.object()["coordinate_provider"]="causal-history";
        meta.object()["value_space"]="history-coordinates";
        meta.object()["coordinate_transform_identity"]=transform.identity();
        meta.object()["coordinate_transform_archive"]=std::filesystem::absolute(path).generic_u8string();
        meta.object()["history_span"]=double(transform.basis().delays-1)*transform.sampleInterval();
        meta.object()["derivative_evidence"]=modelRates?"frozen-projection-of-analytic-measurement-rate":"none";
    } else if(kind!="identity"&&kind!="observable")throw std::invalid_argument("Unknown coordinate transform kind");
    const auto& z=result.coordinates;const auto n=z.sampleCount(),d=z.componentCount();
    c::NpzRecord archive;
    archive.arrays["time"]=c::NpyArray::fromDoubles({n},{z.timeData(),z.timeData()+n});
    archive.arrays["values"]=c::NpyArray::fromDoubles({n,d},{z.valuesData(),z.valuesData()+n*d});
    if(result.derivatives)archive.arrays["coordinate_derivatives"]=c::NpyArray::fromDoubles({n,d},
        {result.derivatives->valuesData(),result.derivatives->valuesData()+n*d});
    meta.object()["coordinate_transform_kind"]=kind;
    meta.object()["coordinate_dimension"]=double(d);
    meta.object()["delta_t"]=n>1&&z.uniformTimeSteps()?z.deltaTime():0;
    archive.arrays["metadata_json"]=c::NpyArray::fromString(json::dump(meta));
    return publish(root,source,"coordinate-transform",label,std::move(archive));
}
}
