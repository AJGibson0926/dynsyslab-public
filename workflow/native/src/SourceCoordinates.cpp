#include "PostprocessInternal.hpp"
#include <cmath>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
bool storedInDimensionalUnits(const json::Value& meta) {
    if(meta.at("value_space").string()=="history-coordinates"||
        meta.at("coordinate_provider").string()=="measured")return false;
    const auto identifier=meta.at("system").string();
    if(identifier!="rpe"&&identifier!="rayleigh-plesset"&&identifier!="marmottant")return false;
    const auto stored=meta.at("archive_coordinate_system").string(
        meta.at("stored_coordinate_system").string());
    const bool legacy=meta.at("execution_backend").at("name").string()=="python-reference"
        ||meta.at("time_series_format").string()=="dynsyslab.TimeSeries/v1";
    return stored=="dimensional"||(stored.empty()&&legacy);
}

void normalizeSourceCoordinates(SourceData& source) {
    if(source.metadataJson.empty())return;
    auto meta=json::parse(source.metadataJson);
    if(!storedInDimensionalUnits(meta))return;
    if(source.dimension!=2)throw std::invalid_argument("bubble archive conversion requires two state components");
    const auto configuration=systemConfiguration(source);
    const auto system=dynsyslab::core::createSystem(configuration);
    const double radius=system->configuration().parameters.at("equilibrium_radius_m");
    const double timeScale=meta.at("nondimensionalization").at("time_scale_s").number(
        meta.at("time_scale_to_seconds").number(system->coordinates().secondsPerModelTime));
    if(!(radius>0)||!(timeScale>0)||!std::isfinite(timeScale))
        throw std::invalid_argument("bubble archive lacks valid coordinate scales");
    for(auto& t:source.time)t/=timeScale;
    for(std::size_t row=0;row<source.samples;++row) {
        source.values[2*row]/=radius;
        source.values[2*row+1]*=timeScale/radius;
    }
    if(!source.inputs.empty()) {
        const auto& history=meta.at("control_history");
        const auto convention=history.at("coordinate_system").string();
        const auto units=history.at("units").string();
        if(convention=="dimensional"||units=="Pa") {
            const auto ambient=system->configuration().parameters.at("ambient_pressure_pa");
            for(auto& input:source.inputs)input/=ambient;
        } else if(convention!="dimensionless"&&units!="p_control/p0")
            throw std::invalid_argument("legacy controlled archive must declare its input units");
    }
    meta.at("system")=system->descriptor().identifier;
    meta.at("source_archive_coordinate_system")="dimensional";
    meta.at("archive_coordinate_system")="dimensionless";
    meta.at("integration_form")="dimensionless";
    meta.at("source_coordinate_transform")=json::Value::Object{
        {"schema","dynsyslab.coordinate-transform/v1"},{"time_multiplier",1/timeScale},
        {"state_multipliers",json::Value::Array{1/radius,timeScale/radius}}};
    json::Value::Object resolved;
    for(const auto& item:system->configuration().parameters)resolved[item.first]=item.second;
    meta.at("resolved_system_parameters")=std::move(resolved);
    meta.at("system_preset")=system->configuration().preset;
    source.metadataJson=json::dump(meta);
}
}  // namespace dynsyslab::workflow::postprocess_detail
