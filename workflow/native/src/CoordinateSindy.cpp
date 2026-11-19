#include "PostprocessInternal.hpp"
#include "KromDatasetSupport.hpp"
#include "ObservationSource.hpp"
#include "DerivativeSettings.hpp"
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
        if(d.shape.size()!=2||d.shape[1]!=source.dimension)
            throw std::invalid_argument("Measured coordinate derivatives must have one column per coordinate");
        if(d.shape[0]==source.samples)derivatives=d.doubles();
        else {  // rows selected by sample_stride / time_window: match the kept samples by time
            const auto all=d.doubles();const auto fullTime=a.at("time").doubles();
            if(fullTime.size()!=d.shape[0])throw std::invalid_argument("Measured coordinate derivatives must align with the archive's time");
            std::size_t j=0;
            for(std::size_t i=0;i<source.samples;++i) {
                while(j<fullTime.size()&&fullTime[j]<source.time[i]-1e-12*std::max(1.0,std::abs(source.time[i])))++j;
                if(j==fullTime.size()||std::abs(fullTime[j]-source.time[i])>1e-9*std::max(1.0,std::abs(source.time[i])))
                    throw std::invalid_argument("Measured coordinate derivatives do not cover the selected rows");
                derivatives.insert(derivatives.end(),all.begin()+std::ptrdiff_t(j*source.dimension),all.begin()+std::ptrdiff_t((j+1)*source.dimension));
            }
        }
    } else if(method=="finite_difference"&&!settings.at("derivative").isNull()) {
        // Noise-robust coordinate derivatives (see NumericalDifferentiation.hpp) instead of the
        // fit's internal finite differences.
        const auto spec=parseDerivativeSpec(settings.at("derivative"));
        const double dt=uniformDeltaTime(source);
        derivatives.assign(source.samples*source.dimension,0.0);
        for(std::size_t k=0;k<source.dimension;++k) {
            std::vector<double> column(source.samples);
            for(std::size_t i=0;i<source.samples;++i)column[i]=source.values[i*source.dimension+k];
            const auto d=c::differentiate(column,dt,spec).derivative;
            for(std::size_t i=0;i<source.samples;++i)derivatives[i*source.dimension+k]=d[i];
        }
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
    meta.object()["gamma_method"]=method;
    if(!settings.at("derivative").isNull())meta.object()["coordinate_derivative"]=settings.at("derivative");meta.object()["coordinate_provider"]=history?"causal-history":"measured";
    meta.object()["delta_t"]=uniformDeltaTime(source);
    std::map<std::string,c::NpyArray> extra;
    extra["training_time"]=c::NpyArray::fromDoubles({source.samples},source.time);
    extra["training_coordinates"]=c::NpyArray::fromDoubles({source.samples,source.dimension},source.values);
    return publish(root,source,"koopman-with-sindy",label,c::makeKoopmanArchive(model,json::dump(meta),std::move(extra)));
}
}
