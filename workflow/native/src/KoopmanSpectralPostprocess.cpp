#include "HankelAnalysis.hpp"
#include "KoopmanSpectralSettings.hpp"
#include <dynsyslab/core/SpectralKromArchive.hpp>
namespace dynsyslab::workflow::postprocess_detail {
std::filesystem::path runKoopmanSpectral(const SourceData& source,const json::Value& settings,
    const std::filesystem::path& root,const std::string& label,const std::string& method) {
    namespace c=dynsyslab::core;using namespace spectral_settings;
    // Reuse the one shared observable/Hankel sampling contract, including immutable matrix reuse.
    auto sampling=settings;
    if(sampling.at("hankel_archive").isNull()&&sampling.at("observables").isNull()&&sampling.at("observable").isNull()) {
        json::Value::Array components;for(std::size_t i=0;i<source.dimension;++i)components.emplace_back(double(i));
        sampling.at("observables")=components;
    }
    if(sampling.at("hankel_archive").isNull()&&sampling.at("window").isNull())
        sampling.at("window")=json::Value::Object{{"delays",1.}};
    auto input=prepareHankel(source,sampling,true);
    const auto& h=input.matrix;c::KoopmanSnapshotPairs data;
    data.samples=h.columns;data.dimension=h.rows;data.x.resize(h.x.size());data.y.resize(h.y.size());
    for(std::size_t i=0;i<data.samples;++i)for(std::size_t j=0;j<data.dimension;++j) {
        data.x[i*data.dimension+j]=h.x[j*data.samples+i];data.y[i*data.dimension+j]=h.y[j*data.samples+i];
    }
    data.weights=numbers(settings.at("weights"));
    if(method=="specrkhs") {
        c::SpectralData result;result.rkhsAdjoints.push_back(c::rkhsAdjointSpectrum(data,settings.at("bandwidth").number(1),
            count(settings.at("rank")),settings.at("svd_tolerance").number(1e-10),complexes(settings.at("residual_grid"))));
        c::NpzRecord archive;archive.arrays["spectral_data_json"]=c::NpyArray::fromString(json::dump(c::spectralDataJson(result)));
        archive.arrays["metadata_json"]=c::NpyArray::fromString(metadata(source,"specrkhs",label,""));
        archiveHankelLayout(archive,input,settings.at("save_matrix").boolean());
        return publish(root,source,"koopman-specrkhs",label,std::move(archive));
    }
    data.successorReplicates=numbers(settings.at("successor_replicates"));data.replicateCount=count(settings.at("replicate_count"));
    auto fitSettings=settings;if(method!="koopman_spectral")fitSettings.at("method")=method;
    auto fitOptions=options(fitSettings,uniformDeltaTime(source)*h.window.sampleStride);
    if(!settings.at("sample_interval").isNull()&&std::abs(fitOptions.sampleInterval-uniformDeltaTime(source)*h.window.sampleStride)>1e-10*fitOptions.sampleInterval)
        throw std::invalid_argument("Spectral sample interval differs from paired Hankel shift");
    auto model=c::fitKoopmanSpectral(data,dictionary(settings,data),fitOptions);
    model.coordinateContract=json::dump(input.provenance);
    c::fitSpectralReconstruction(model,"coordinates",data.x,data.x,data.dimension,data.weights);
    measures(model,settings);
    auto archive=c::makeSpectralKromArchive(model,metadata(source,method,label,""));
    archiveHankelLayout(archive,input,settings.at("save_matrix").boolean());
    archive.arrays["eigenfunction_values"]=c::NpyArray::fromComplex({data.samples,model.rank},model.evaluateEigenfunctions(data.x));
    return publish(root,source,"koopman-"+method,label,std::move(archive));
}
}
