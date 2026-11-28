// The koopman_spectral family of post-processing methods (EDMD, ResDMD, mpEDMD, Hermitian, TLS and
// multiplicative DMD, stochastic ResDMD; Williams, Kevrekidis & Rowley 2015; Colbrook & Townsend 2024;
// Colbrook 2023) and SpecRKHS. A fit is published as one krom: its eigenpairs with their residuals, the
// pseudospectrum, its compressed operator, and two representations of the same eigenfunctions, the
// frozen-dictionary expansion and an interpolation over the physical states at the training columns.
// Settings "additional_sources" pool further records of the same system into the fit: each record's
// snapshot pairs carry its transverse weight (settings "source_weights", OrbitAverages.hpp: uniform, a
// list, or canonical area) divided by its number of pairs, so a family of orbits discretizes the
// invariant measure (time averages along each orbit, the weights across them).
#include "HankelAnalysis.hpp"
#include "TimeSeriesSource.hpp"
#include "OrbitAverages.hpp"
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <algorithm>
#include <cmath>
#include <complex>
#include "KoopmanSpectralSettings.hpp"
namespace dynsyslab::workflow::postprocess_detail {
namespace {
namespace c=dynsyslab::core;
// The domain of the frozen-dictionary eigenfunctions, from the Hankel window: one delay of a solve's
// full state is the physical state; one delay of named coordinate columns is those catalog observables;
// longer windows are data coordinates.
c::FunctionDomain spectralDomain(const SourceData& source,const json::Value& window,std::size_t dimension,
                                 const std::string& method,const std::string& label) {
    c::FunctionDomain domain;domain.dimension=dimension;
    const auto& channels=window.at("channels");
    bool oneDelay=window.at("delays").number(0)==1&&!channels.isNull()&&channels.array().size()==dimension;
    std::vector<std::string> observables;
    for(std::size_t k=0;oneDelay&&k<channels.array().size();++k) {
        const auto& channel=channels.array()[k];
        oneDelay=channel.at("lag_samples").number(0)==0&&channel.at("measurement_column").number(-1)==double(k);
        observables.push_back(channel.at("observable").string());
    }
    if(!oneDelay)return domain;
    const auto meta=json::parse(metadata(source,method,label,""));
    const bool fromSolve=meta.at("coordinate_transform_kind").string().empty()&&meta.at("value_space").string().empty();
    std::size_t physical=0;
    try{physical=c::createSystem(systemConfiguration(source))->descriptor().stateDimension;}catch(const std::exception&){}
    if(fromSolve&&physical==dimension)return c::FunctionDomain::physicalState(dimension);
    if(!fromSolve&&std::none_of(observables.begin(),observables.end(),[](const std::string& n){return n.empty();}))
        return c::FunctionDomain::fromCatalog(observables,dimension);
    return domain;
}
}
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
    // The lead record and any pooled records, each with its Hankel pairs.
    std::vector<SourceData> pooled{source};
    for(const auto& path:settings.at("additional_sources").array())pooled.push_back(loadSource(std::filesystem::u8path(path.string())));
    std::vector<HankelAnalysis> inputs{input};
    for(std::size_t s=1;s<pooled.size();++s) {
        inputs.push_back(prepareHankel(pooled[s],sampling,true));
        if(inputs.back().matrix.rows!=h.rows||json::dump(inputs.back().provenance.at("channels"))!=json::dump(input.provenance.at("channels"))
           ||std::abs(uniformDeltaTime(pooled[s])-uniformDeltaTime(source))>1e-12*uniformDeltaTime(source))
            throw std::invalid_argument("pooled records must give the same delay coordinates on the same clock");
    }
    data.dimension=h.rows;
    for(const auto& in:inputs) {
        const auto& m=in.matrix;
        for(std::size_t i=0;i<m.columns;++i)for(std::size_t j=0;j<m.rows;++j) {
            data.x.push_back(m.x[j*m.columns+i]);data.y.push_back(m.y[j*m.columns+i]);
        }
        data.samples+=m.columns;
    }
    data.weights=numbers(settings.at("weights"));
    if(pooled.size()>1) {
        if(!data.weights.empty())throw std::invalid_argument("pooled fits weight records with source_weights, not per-pair weights");
        const auto w=recordWeights(settings.at("source_weights"),pooled);
        for(std::size_t s=0;s<inputs.size();++s)
            for(std::size_t i=0;i<inputs[s].matrix.columns;++i)data.weights.push_back(w[s]/double(inputs[s].matrix.columns));
    }
    c::NpzRecord layout;archiveHankelLayout(layout,input,settings.at("save_matrix").boolean());
    if(method=="specrkhs") {  // adjoint sections, not a forward model: a report
        c::SpectralDiagnostics result;result.rkhsAdjoints.push_back(c::rkhsAdjointSpectrum(data,settings.at("bandwidth").number(1),
            count(settings.at("rank")),settings.at("svd_tolerance").number(1e-10),complexes(settings.at("residual_grid"))));
        layout.arrays["spectral_diagnostics_json"]=c::NpyArray::fromString(json::dump(c::spectralDiagnosticsJson(result)));
        layout.arrays["metadata_json"]=c::NpyArray::fromString(metadata(source,"specrkhs",label,""));
        return publishObject(root, source, dynsyslab::core::object_kind::report, "koopman-specrkhs",std::move(layout));
    }
    if(!settings.at("rigged_measures").isNull())
        throw std::invalid_argument("rigged_measures is a stage of its own: add a rigged_measure stage on this krom");
    data.successorReplicates=numbers(settings.at("successor_replicates"));data.replicateCount=count(settings.at("replicate_count"));
    auto fitSettings=settings;if(method!="koopman_spectral")fitSettings.at("method")=method;
    auto fitOptions=options(fitSettings,uniformDeltaTime(source)*h.window.sampleStride);
    if(!settings.at("sample_interval").isNull()&&std::abs(fitOptions.sampleInterval-uniformDeltaTime(source)*h.window.sampleStride)>1e-10*fitOptions.sampleInterval)
        throw std::invalid_argument("Spectral sample interval differs from paired Hankel shift");
    auto fit=c::fitKoopmanSpectral(data,dictionary(settings,data),fitOptions,
                                   spectralDomain(source,input.provenance,data.dimension,method,label));
    {
        auto compression=*fit.compression();compression.coordinateContract=json::dump(input.provenance);
        fit.setCompression(std::move(compression));
    }
    c::fitFeatureReadout(fit,"coordinates",data.x,data.x,data.dimension,data.weights);
    // The eigenfunctions' values at the training columns and the physical states there, for the
    // interpolation representation.
    const double dt=uniformDeltaTime(source)*double(h.window.sampleStride);
    const auto multipliers=fit.spectrum().discreteMultipliers();
    std::vector<std::complex<double>> rates;
    for(const auto mu:multipliers)rates.push_back(std::log(mu)/dt);
    const auto& expansion=c::frozenDictionaryExpansion(fit);
    c::EigenfunctionSamples samples;
    samples.modes=multipliers.size();
    {
        const auto phi=fit.compression()->dictionary.evaluate(data.x);const auto n=fit.compression()->dictionary.size();
        samples.values.assign(data.samples*samples.modes,{});
        for(std::size_t i=0;i<data.samples;++i)for(std::size_t j=0;j<samples.modes;++j)for(std::size_t k=0;k<n;++k)
            samples.values[i*samples.modes+j]+=phi[i*n+k]*expansion.coefficients()[j*n+k];
    }
    // The physical states at every training column, record by record (sample indices count across the
    // pooled records in order).
    std::vector<double> states;
    std::size_t width=0,offset=0;
    for(std::size_t s=0;s<pooled.size();++s) {
        const auto record=c::readNpz(pooled[s].path);
        const bool carried=record.contains("source_states");
        const auto physical=carried?record.at("source_states").doubles():pooled[s].values();
        width=carried?record.at("source_states").shape.at(1):pooled[s].dimension();
        const auto& m=inputs[s].matrix;
        for(std::size_t j=0;j<m.columns;++j) {
            const auto tick=inputs[s].alignEnd?m.ends.at(j):m.starts.at(j);
            samples.time.push_back(pooled[s].time().at(tick));
            samples.sampleIndices.push_back(offset+tick);
            states.insert(states.end(),physical.begin()+std::ptrdiff_t(tick*width),physical.begin()+std::ptrdiff_t((tick+1)*width));
        }
        offset+=pooled[s].samples();
    }
    std::vector<std::shared_ptr<const c::FunctionRepresentation>> representations;
    representations.push_back(std::make_shared<const c::FunctionExpansion>(expansion));
    representations.push_back(c::interpolateSamples(states,data.samples,c::FunctionDomain::physicalState(width),samples.values,{},8,
                                                    method+" eigenfunctions at their training columns"));
    c::DataOrigin origin{method,settings,{dataReference(source.path,"time_series")}};
    for(std::size_t s=1;s<pooled.size();++s)origin.inputs.push_back(dataReference(pooled[s].path,"time_series"));
    // The K-ROM evaluates the first representation: the frozen expansion when it is a function of the
    // state (one delay of a solve's state), else the interpolation on the physical states.
    if(representations.front()->coordinates().requiresExternalValues())
        std::swap(representations.front(),representations.back());
    auto spectrum=c::Spectrum::fromDiscrete(multipliers,rates,dt,origin);
    // Residuals, the pseudospectrum with its minimizers (Colbrook & Townsend 2024), measures and packets.
    spectrum.adoptEvidence(fit.spectrum());
    c::KoopmanModel krom(std::move(spectrum),std::move(samples),origin,std::move(representations));
    krom.setCompression(*fit.compression());
    krom.setSpectralDiagnostics(fit.spectralDiagnostics());  // invariance, stochastic residuals
    layout.arrays.erase("metadata_json");
    return publishKrom(root,source,krom,"koopman-"+method,metadata(source,method,label,""),std::move(layout.arrays),source.record.channels);
}
}
