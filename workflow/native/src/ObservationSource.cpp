#include "ObservationSource.hpp"
#include <dynsyslab/core/ScientificCompatibility.hpp>
#include <dynsyslab/core/TimeSeries.hpp>
#include <cmath>
#include <stdexcept>
namespace dynsyslab::workflow::postprocess_detail {
SourceData loadObservationSource(const std::filesystem::path& path,
    const std::string& arrayName,const std::string& analyticObservable) {
    const auto a=dynsyslab::core::readNpz(path);
    if(!a.contains(arrayName)) {
        if(analyticObservable.empty())throw std::invalid_argument("Missing measured array: "+arrayName);
        auto s=loadSource(path);
        const auto system=dynsyslab::core::createSystem(systemConfiguration(s));
        const auto observable=system->observableCatalog().at(analyticObservable);
        return withSeries(s,s.time(),observable.evaluateRealBatch(s.values(),s.samples()),observable.outputDimension());
    }
    const auto& t=a.at("time");const auto& y=a.at(arrayName);
    if(t.shape.size()!=1||(y.shape.size()!=1&&y.shape.size()!=2)||y.shape[0]!=t.shape[0]||t.shape[0]<2)
        throw std::invalid_argument("Measured data require a samples-by-components array aligned with time");
    SourceData s{path,a.metadataJson(),{}};
    s.record.series=dynsyslab::core::TimeSeries(t.doubles(),y.doubles(),y.shape.size()==1?1:y.shape[1]);
    for(double yv:s.values())if(!std::isfinite(yv))throw std::invalid_argument("Measured values must be finite");
    const std::string key=a.contains("control_interval_values")?"control_interval_values":"control_values";
    if(a.contains(key)) {
        const auto& u=a.at(key);
        if(u.shape.size()!=2||!u.shape[1]||u.shape[0]!=s.samples()-(key=="control_interval_values"?1:0))
            throw std::invalid_argument("Measured input history has incompatible shape");
        s.record.inputsOnIntervals=key=="control_interval_values";
        s.record.inputs=inputSeries(s.time(),u.doubles(),u.shape[1],s.record.inputsOnIntervals);
        for(double v:s.inputs())if(!std::isfinite(v))throw std::invalid_argument("Measured inputs must be finite");
    }
    return s;
}
std::optional<dynsyslab::core::SystemConfiguration> observationSystem(const SourceData& s) {
    const auto meta=s.metadataJson.empty()?json::Value{}:json::parse(s.metadataJson);
    if(meta.at("system").string().empty()&&meta.at("system_configuration").at("identifier").string().empty())return {};
    return systemConfiguration(s);
}
void requireSameObservationSystem(const SourceData& expected,const SourceData& actual) {
    const auto e=observationSystem(expected),a=observationSystem(actual);
    if(bool(e)!=bool(a))throw std::invalid_argument("Cannot verify a pool mixing declared and undeclared physical systems");
    if(e)dynsyslab::core::requireCompatibleSystems(*e,*a);
    const auto em=expected.metadataJson.empty()?json::Value{}:json::parse(expected.metadataJson);
    const auto am=actual.metadataJson.empty()?json::Value{}:json::parse(actual.metadataJson);
    for(const char* key:{"measurement_time_units","measurement_units"})
        if(em.at(key).string()!=am.at(key).string())throw std::invalid_argument(std::string("Measurement unit mismatch: ")+key);
}
}
