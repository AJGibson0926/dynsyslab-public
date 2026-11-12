#include "PostprocessInternal.hpp"
#include "KromDatasetSupport.hpp"
#include "ObservationSource.hpp"
#include <dynsyslab/core/EigenfunctionFit.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace dynsyslab::workflow::postprocess_detail {
std::filesystem::path runCoordinateSindy(const SourceData& source,const json::Value& settings,
    const std::filesystem::path& spectrumPath,const std::vector<std::complex<double>>& rates,
    const dynsyslab::core::KoopmanSINDyOptions& options,
    const std::filesystem::path& root,const std::string& label) {
    namespace c=dynsyslab::core;
    if(datasetDynamics(source)!="unforced")throw std::invalid_argument("Autonomous coordinate SINDy requires confirmed unforced training data");
    const auto metadataSource=json::parse(source.metadataJson);
    c::EigenfunctionCoordinates q;q.dimension=source.dimension;
    const bool history=metadataSource.at("coordinate_provider").string()=="causal-history";
    q.kind=history?c::EigenfunctionCoordinateKind::causalHistory:c::EigenfunctionCoordinateKind::measured;
    q.providerIdentity=metadataSource.at("coordinate_transform_identity").string();
    q.providerArtifact=metadataSource.at("coordinate_transform_archive").string();
    if(history&&q.providerIdentity.empty())throw std::invalid_argument("History-coordinate fitting needs the frozen transform identity");
    for(std::size_t k=0;k<source.dimension;++k)q.names.push_back((history?"xi_":"q_")+std::to_string(k+1));
    const auto method=settings.at("gamma_method").string("finite_difference");
    std::vector<double> derivatives;
    if(method=="measured") {
        const auto a=c::readNpz(source.path);
        const auto& d=a.at(settings.at("derivative_array").string("coordinate_derivatives"));
        if(d.shape!=std::vector<std::size_t>{source.samples,source.dimension})
            throw std::invalid_argument("Measured coordinate derivatives must align with the snapshots");
        derivatives=d.doubles();
    } else if(method!="finite_difference")throw std::invalid_argument("Measured/history SINDy uses measured derivatives or finite differences, never a physical RHS");
    auto fit=c::fitSindyEigenfunctions(source.values,source.samples,source.time,rates,q,derivatives,options);
    auto physical=observationSystem(source);
    const double inputs=settings.at("input_dimension").number(double(source.inputDimension));
    if(!std::isfinite(inputs)||inputs<0||inputs!=std::floor(inputs)
        ||inputs>double(std::numeric_limits<std::size_t>::max()/2))
        throw std::invalid_argument("Coordinate SINDy input_dimension must be a nonnegative integer");
    auto model=physical?c::KoopmanModel(fit.eigenfunctions,*physical,fit.diagnostics):
        c::KoopmanModel(fit.eigenfunctions,static_cast<std::size_t>(inputs),fit.diagnostics);
    if(physical&&!settings.at("input_dimension").isNull()&&inputs!=double(model.inputDimension()))
        throw std::invalid_argument("Declared input dimension conflicts with the physical system");
    model.addDatasetReference({std::filesystem::absolute(source.path).generic_u8string(),"training","unforced",{}});
    auto meta=json::parse(metadata(source,"koopman-with-sindy",label,""));
    meta.object()["source_spectrum_archive"]=spectrumPath.generic_u8string();
    meta.object()["gamma_method"]=method;meta.object()["coordinate_provider"]=history?"causal-history":"measured";
    meta.object()["delta_t"]=uniformDeltaTime(source);
    std::map<std::string,c::NpyArray> extra;
    extra["training_time"]=c::NpyArray::fromDoubles({source.samples},source.time);
    extra["training_coordinates"]=c::NpyArray::fromDoubles({source.samples,source.dimension},source.values);
    return publish(root,source,"koopman-with-sindy",label,c::makeKoopmanArchive(model,json::dump(meta),std::move(extra)));
}
}
