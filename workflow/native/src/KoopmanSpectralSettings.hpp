#pragma once
#include <dynsyslab/core/KoopmanSpectralAnalysis.hpp>
#include <dynsyslab/core/Archive.hpp>
#include <cmath>
#include <stdexcept>
#include <limits>
namespace dynsyslab::workflow::spectral_settings {
namespace c=dynsyslab::core;using J=c::json::Value;
inline std::size_t count(const J& v,std::size_t fallback=0) {
    if(v.isNull())return fallback;const double n=v.number(-1);
    if(!std::isfinite(n)||n<0||n>9007199254740991.||std::floor(n)!=n)throw std::invalid_argument("Expected nonnegative integer count");
    return std::size_t(n);
}
inline std::vector<double> numbers(const J& a) {
    std::vector<double> out;
    for(const auto& v:a.array()) {
        const double n=v.number(std::numeric_limits<double>::quiet_NaN());
        if(!std::isfinite(n))throw std::invalid_argument("Expected finite flat numeric array");out.push_back(n);
    }
    return out;
}
inline std::vector<std::complex<double>> complexes(const J& a) {
    std::vector<std::complex<double>> out;for(const auto& v:a.array()) {
        auto z=numbers(v);if(z.size()!=2)throw std::invalid_argument("Complex values require [real,imag]");out.emplace_back(z[0],z[1]);
    }return out;
}
// Pseudospectrum points: a list of [real, imag], or a rectangular lattice {"real": [from, to, n],
// "imag": [from, to, n]} (endpoints included; real varies fastest), which the Viewer draws as a map.
inline std::vector<std::complex<double>> residualGrid(const J& v) {
    if(!std::holds_alternative<J::Object>(v.storage()))return complexes(v);
    const auto re=numbers(v.at("real")),im=numbers(v.at("imag"));
    if(re.size()!=3||im.size()!=3||!(re[2]>=2)||!(im[2]>=2)||re[2]*im[2]>4e6)
        throw std::invalid_argument("A residual lattice needs real and imag as [from, to, n], n >= 2");
    std::vector<std::complex<double>> out;
    const auto nr=std::size_t(re[2]),ni=std::size_t(im[2]);
    for(std::size_t j=0;j<ni;++j)for(std::size_t k=0;k<nr;++k)
        out.emplace_back(re[0]+(re[1]-re[0])*double(k)/double(nr-1),im[0]+(im[1]-im[0])*double(j)/double(ni-1));
    return out;
}
// Spectral angles: a list, or {"count": n} for n uniform angles on [-pi, pi) (optionally "from", "to",
// endpoint excluded, and "append": [extra angles]).
inline std::vector<double> angleGrid(const J& v) {
    if(std::holds_alternative<J::Array>(v.storage()))return numbers(v);
    const double pi=3.14159265358979323846;const auto n=count(v.at("count"));
    const double from=v.at("from").number(-pi),to=v.at("to").number(pi);
    if(!n||n>1000000||!std::isfinite(from)||!std::isfinite(to)||!(to>from))throw std::invalid_argument("An angle grid needs a count and from < to");
    std::vector<double> out;for(std::size_t k=0;k<n;++k)out.push_back(from+(to-from)*double(k)/double(n));
    if(!v.at("append").isNull())for(double x:numbers(v.at("append")))out.push_back(x);
    return out;
}
inline J array(const std::vector<double>& a){J::Array out;for(double v:a)out.emplace_back(v);return out;}
inline J array(const std::vector<std::complex<double>>& a){J::Array out;for(auto v:a)out.emplace_back(J::Array{v.real(),v.imag()});return out;}
inline c::KoopmanSpectralOptions options(const J& settings,double dt=1) {
    c::KoopmanSpectralOptions o;o.method=settings.at("method").string("resdmd");o.rank=count(settings.at("rank"));
    o.relativeTolerance=settings.at("svd_tolerance").number(1e-10);o.sampleInterval=settings.at("sample_interval").number(dt);
    o.acknowledgeMeasurePreserving=settings.at("assume_measure_preserving").boolean();
    o.acknowledgeSelfAdjoint=settings.at("assume_self_adjoint").boolean();
    o.weightPolicy=settings.at("weight_policy").string("uniform");
    o.residualGrid=residualGrid(settings.at("residual_grid"));return o;
}
inline c::SpectralDictionary dictionary(const J& settings,const c::KoopmanSnapshotPairs& data) {
    auto d=settings.at("dictionary");
    if(!d.isNull()&&!std::holds_alternative<J::Object>(d.storage()))throw std::invalid_argument("Spectral dictionary must be an object with kind and settings");
    const auto& requested=static_cast<const J&>(d);
    if(!requested.at("nodes").isNull()||!requested.at("exponents").isNull()) {
        d.at("dimension")=double(data.dimension);return c::spectralDictionaryFromJson(d);
    }
    c::SpectralDictionaryOptions o;o.kind=requested.at("kind").string("linear");
    o.degree=count(requested.at("degree"),2);o.centers=count(requested.at("centers"),32);
    o.maximumFeatures=count(requested.at("maximum_features"),4096);o.constant=requested.at("constant").boolean(true);
    o.normalize=requested.at("normalize").boolean();o.bandwidth=requested.at("bandwidth").number(1);
    o.indexSet=requested.at("index_set").string("total_degree");
    o.centerSelection=requested.at("center_selection").string("farthest");
    o.centerIterations=count(requested.at("center_iterations"),100);
    o.angularOrder=count(requested.at("angular_order"),4);o.hermiteOrder=count(requested.at("hermite_order"),4);
    o.smoothness=requested.at("smoothness").number(2);o.amplitude=requested.at("amplitude").number(1);
    o.center=numbers(requested.at("center"));o.scale=numbers(requested.at("scale"));
    return c::buildSpectralDictionary(data.x,data.dimension,o);
}
inline c::KoopmanSnapshotPairs pairs(const J& value) {
    if(!value.at("archive").isNull()) {
        if(!value.at("x").isNull()||!value.at("y").isNull())throw std::invalid_argument("Specify inline pairs or an archive, not both");
        const auto record=c::readNpz(std::filesystem::u8path(value.at("archive").string()));
        const auto metadata=c::json::parse(record.metadataJson());
        if(metadata.at("schema").string()!="dynsyslab.koopman-snapshot-pairs/v1")throw std::invalid_argument("Expected versioned Koopman snapshot-pairs archive");
        const auto& x=record.at("x");const auto& y=record.at("y");
        if(x.shape.size()!=2||x.shape!=y.shape||!x.shape[0]||!x.shape[1])throw std::invalid_argument("Snapshot archive x/y need equal nonempty sample-by-coordinate shapes");
        c::KoopmanSnapshotPairs data;data.samples=x.shape[0];data.dimension=x.shape[1];data.x=x.doubles();data.y=y.doubles();
        if(record.arrays.count("weights")){
            const auto& w=record.at("weights");if(w.shape!=std::vector<std::size_t>{data.samples})throw std::invalid_argument("Snapshot weights shape mismatch");data.weights=w.doubles();
        }
        return data;
    }
    c::KoopmanSnapshotPairs data;data.dimension=count(value.at("dimension"));data.x=numbers(value.at("x"));data.y=numbers(value.at("y"));
    if(!data.dimension||data.x.size()%data.dimension)throw std::invalid_argument("Invalid snapshot coordinates");
    data.samples=data.x.size()/data.dimension;data.weights=numbers(value.at("weights"));
    data.successorReplicates=numbers(value.at("successor_replicates"));data.replicateCount=count(value.at("replicate_count"));return data;
}
// Rigged DMD of each listed observable on a unitary compression (Colbrook, Drysdale & Horning 2024):
// [{name, angles, bandwidth, order, coefficients | reconstruction{name, output}}]. Each appends a
// SpectralMeasure, with its generalized eigenfunctions, to model.spectral.measures.
// Rigged DMD of each listed observable. An item names its band ("band", default "name") and
// its observable: raw dictionary "coefficients" (observable "observable", default "name"), or output
// "output" of a fitted feature "reconstruction" {name, output} (observable: the readout's name, whose
// Koopman modes are set for all its outputs when it has none).
inline void riggedMeasures(c::KoopmanModel& model,const J& list) {
    if(!model.compression())throw std::invalid_argument("Rigged DMD needs a krom with a compressed operator (a spectral fit's)");
    for(const auto& v:list.array()) {
        const auto name=v.at("name").string();
        const auto band=v.at("band").string(name);
        auto coefficients=numbers(v.at("coefficients"));std::string projection,observable=v.at("observable").string(name);
        std::size_t output=0;
        if(!v.at("reconstruction").isNull()) {
            if(!v.at("coefficients").isNull())throw std::invalid_argument("Choose raw coefficients or a fitted reconstruction, not both");
            const auto& ref=v.at("reconstruction");observable=ref.at("name").string();
            const auto& readout=model.compression()->readouts.at(observable);output=count(ref.at("output"));
            if(output>=readout.outputs)throw std::invalid_argument("Reconstruction output index out of range");
            if(!model.spectrum().koopmanModes(observable))
                c::setDictionaryObservableModes(model,observable,readout.coefficients,readout.outputs);
            coefficients.resize(model.compression()->dictionary.size());
            for(std::size_t i=0;i<coefficients.size();++i)coefficients[i]=readout.coefficients[i*readout.outputs+output];
            projection="Observable explicitly selected from fitted reconstruction '"+observable+"', output "+std::to_string(output)+"; the readout's relativeTrainingError is for all its outputs. This measures the retained projection, which can differ from the original observable.";
        }
        c::riggedSpectralMeasure(model,coefficients,observable,output,band,
            angleGrid(v.at("angles")),v.at("bandwidth").number(),count(v.at("order"),1));
        if(!projection.empty()) {  // record the explicit projection on the measure just added
            auto spectrum=model.spectrum();auto measures=spectrum.measures();
            for(auto& m:measures)if(m.band==band&&m.observable==observable&&m.output==output)m.assumptions=projection;
            spectrum.setMeasures(std::move(measures));
            model.setSpectralEvidence(spectrum);
        }
    }
}
} // namespace dynsyslab::workflow::spectral_settings
