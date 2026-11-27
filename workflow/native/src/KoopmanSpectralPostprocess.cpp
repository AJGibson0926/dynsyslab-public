#include "HankelAnalysis.hpp"
#include "TimeSeriesSource.hpp"
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <cmath>
#include <complex>
#include "KoopmanSpectralSettings.hpp"
#include <dynsyslab/core/SpectralKromArchive.hpp>
namespace dynsyslab::workflow::postprocess_detail {
std::filesystem::path runKoopmanSpectral(const SourceData& source,const json::Value& settings,
    const std::filesystem::path& root,const std::string& label,const std::string& method) {
    namespace c=dynsyslab::core;using namespace spectral_settings;
    // Reuse the one shared observable/Hankel sampling contract, including immutable matrix reuse.
    auto sampling=settings;
    if(sampling.at("hankel_archive").isNull()&&sampling.at("observables").isNull()&&sampling.at("observable").isNull()) {
        json::Value::Array components;for(std::size_t i=0;i<source.dimension();++i)components.emplace_back(double(i));
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
        return publishObject(root, source, dynsyslab::core::object_kind::report, "koopman-specrkhs",std::move(archive));
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
    // The fit's SpectralData is kept once, on the krom: residuals on its eigenpairs, the residual grid as
    // its pseudospectrum, measures as continuous components (Spectrum.hpp), the rest as its spectral
    // diagnostics. The embedded operator contract carries none of it.
    const auto evidence=std::move(model.spectral);
    model.spectral={};
    auto archive=c::makeSpectralKromArchive(model,metadata(source,method,label,""));
    archiveHankelLayout(archive,input,settings.at("save_matrix").boolean());
    // The result is a krom (KoopmanArchive.hpp) with two representations of the same eigenfunctions:
    // the frozen-dictionary expansion, and an interpolation over the physical states at the training
    // columns. The pairs' residuals describe the frozen-dictionary expansion the method fitted.
    const double dt=uniformDeltaTime(source)*double(h.window.sampleStride);
    std::vector<std::complex<double>> rates;
    for(const auto mu:model.multipliers)rates.push_back(std::log(mu)/dt);
    const auto record=c::readNpz(source.path);
    const bool carried=record.contains("source_states");
    const auto physical=carried?record.at("source_states").doubles():source.values();
    const std::size_t width=carried?record.at("source_states").shape.at(1):source.dimension();
    c::EigenfunctionSamples samples;
    samples.modes=model.rank;
    samples.values=model.evaluateEigenfunctions(data.x);
    std::vector<double> states;
    for(std::size_t j=0;j<h.columns;++j) {
        const auto tick=input.alignEnd?h.ends.at(j):h.starts.at(j);
        samples.time.push_back(source.time().at(tick));
        samples.sampleIndices.push_back(tick);
        states.insert(states.end(),physical.begin()+std::ptrdiff_t(tick*width),physical.begin()+std::ptrdiff_t((tick+1)*width));
    }
    std::vector<std::shared_ptr<const c::FunctionRepresentation>> representations;
    try {  // the frozen dictionary, on the domain its window describes (needs the producing system)
        const auto expansion=c::spectralKromModel(archive,systemConfiguration(source));
        representations.push_back(expansion.representation().function);
    } catch(const std::exception&) {}
    representations.push_back(c::interpolateSamples(states,h.columns,c::FunctionDomain::physicalState(width),samples.values,{},8,
                                                    method+" eigenfunctions at their training columns"));
    c::DataOrigin origin{method,settings,{dataReference(source.path,"time_series")}};
    // The K-ROM evaluates the first representation: the frozen expansion when it is a function of the
    // state (one delay of a solve's state), else the interpolation on the physical states.
    if(representations.size()==2&&representations.front()->coordinates().requiresExternalValues())
        std::swap(representations.front(),representations.back());
    auto spectrum=c::Spectrum::fromDiscrete(model.multipliers,rates,dt,origin);
    c::attachSpectralData(spectrum,evidence);  // ResDMD (Colbrook & Townsend 2024): one-step residuals, pseudospectrum
    c::KoopmanModel krom(std::move(spectrum),std::move(samples),origin,std::move(representations));
    krom.setSpectralDiagnostics(evidence);  // invariance, wave packets, stochastic residuals
    const auto meta=archive.metadataJson();
    archive.arrays.erase("metadata_json");
    return publishKrom(root,source,krom,"koopman-"+method,meta,std::move(archive.arrays),source.record.channels);
}
}
