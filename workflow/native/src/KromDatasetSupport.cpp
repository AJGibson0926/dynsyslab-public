#include "KromDatasetSupport.hpp"
#include <dynsyslab/core/ScientificCompatibility.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
std::string datasetDynamics(const SourceData& source) {
    if(std::any_of(source.inputs.begin(),source.inputs.end(),[](double u){return !std::isfinite(u)||u!=0.;}))
        return "forced";
    const auto meta=source.metadataJson.empty()?json::Value{}:json::parse(source.metadataJson);
    const auto forcing=meta.at("forcing_type").string();
    const auto controller=meta.at("controller_kind").string(meta.at("controller_configuration").at("kind").string());
    if((!controller.empty()&&controller!="none") || (!forcing.empty()&&forcing!="free"&&forcing!="unforced"&&forcing!="none"))
        return "forced";
    if(!source.inputs.empty() || controller=="none" || forcing=="free" || forcing=="unforced")
        return "unforced";
    return "unknown";
}

void requireMatchingUnforcedDataset(const dynsyslab::core::KoopmanModel& model,const SourceData& source) {
    if(datasetDynamics(source)!="unforced")
        throw std::invalid_argument("Autonomous Koopman validation requires a confirmed unforced dataset; forced/unknown data cannot be scored with the free semigroup law");
    const auto actual=systemConfiguration(source);
    dynsyslab::core::requireCompatibleSystems(model.systemConfiguration(),actual);
    if(model.requiresExternalCoordinates())
        throw std::invalid_argument("History-coordinate KROM validation requires the frozen coordinate provider; raw dynamical states are not latent coordinates");
    const auto meta=json::parse(source.metadataJson);
    const auto coordinates=meta.at("archive_coordinate_system").string(meta.at("integration_form").string());
    const auto preferred=dynsyslab::core::createSystem(model.systemConfiguration())->descriptor().preferredCoordinates;
    if(!coordinates.empty() && coordinates!="automatic" && coordinates!=preferred)
        throw std::invalid_argument("Validation coordinates differ from the KROM system coordinates");
}

SourceData validationSamples(const SourceData& source,const dynsyslab::core::KoopmanModel& model,const json::Value& settings) {
    const double sourceDt=uniformDeltaTime(source);
    const double requested=settings.at("validation_delta_t").number(
        model.spectrum().sampleInterval()>0?model.spectrum().sampleInterval():sourceDt);
    const double ratio=requested/sourceDt;
    if(!std::isfinite(ratio)||ratio<1-1e-8||ratio>1e9||std::abs(ratio-std::round(ratio))>1e-7*std::max(1.,ratio))
        throw std::invalid_argument("Validation dt must be an integer multiple of the source dt; explicitly resample separately if necessary");
    const auto stride=static_cast<std::size_t>(std::llround(ratio));
    const auto& spectrum=model.spectrum();
    if(spectrum.domain()==dynsyslab::core::EigenvalueDomain::discreteTime) {
        const double multiple=requested/spectrum.sampleInterval();
        if(multiple<1-1e-8 || std::abs(multiple-std::round(multiple))>1e-7)
            throw std::invalid_argument("Sampled-origin eigenfunctions must be validated at their fitted dt or an integer multiple");
    }
    const double start=settings.at("start").number(source.time.front());
    const double stop=settings.at("stop").number(source.time.back());
    if(!std::isfinite(start)||!std::isfinite(stop)||!(stop>start))
        throw std::invalid_argument("Validation interval requires finite start < stop");
    SourceData result=source;result.time.clear();result.values.clear();result.inputs.clear();
    for(std::size_t k=0;k<source.samples;k+=stride) {
        if(source.time[k]<start-1e-10 || source.time[k]>stop+1e-10)continue;
        result.time.push_back(source.time[k]);
        result.values.insert(result.values.end(),source.values.begin()+k*source.dimension,
            source.values.begin()+(k+1)*source.dimension);
    }
    result.samples=result.time.size();result.inputDimension=0;
    if(result.samples<3)throw std::invalid_argument("Validation selection has fewer than three samples");
    return result;
}
}
