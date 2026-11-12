#pragma once
#include <dynsyslab/core/KoopmanSpectralAnalysis.hpp>
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
inline J array(const std::vector<double>& a){J::Array out;for(double v:a)out.emplace_back(v);return out;}
inline J array(const std::vector<std::complex<double>>& a){J::Array out;for(auto v:a)out.emplace_back(J::Array{v.real(),v.imag()});return out;}
inline c::KoopmanSpectralOptions options(const J& settings,double dt=1) {
    c::KoopmanSpectralOptions o;o.method=settings.at("method").string("resdmd");o.rank=count(settings.at("rank"));
    o.relativeTolerance=settings.at("svd_tolerance").number(1e-10);o.sampleInterval=settings.at("sample_interval").number(dt);
    o.acknowledgeMeasurePreserving=settings.at("assume_measure_preserving").boolean();
    o.acknowledgeSelfAdjoint=settings.at("assume_self_adjoint").boolean();
    o.weightPolicy=settings.at("weight_policy").string("uniform");
    o.residualGrid=complexes(settings.at("residual_grid"));return o;
}
inline c::SpectralDictionary dictionary(const J& settings,const c::KoopmanSnapshotPairs& data) {
    auto d=settings.at("dictionary");
    if(!d.isNull()&&!std::holds_alternative<J::Object>(d.storage()))throw std::invalid_argument("Spectral dictionary must be an object with kind and settings");
    const auto& requested=static_cast<const J&>(d);
    if(!requested.at("nodes").isNull()||!requested.at("exponents").isNull()||!requested.at("center").isNull()||!requested.at("scale").isNull()) {
        d.at("dimension")=double(data.dimension);return c::spectralDictionaryFromJson(d);
    }
    c::SpectralDictionaryOptions o;o.kind=requested.at("kind").string("linear");
    o.degree=count(requested.at("degree"),2);o.centers=count(requested.at("centers"),32);
    o.maximumFeatures=count(requested.at("maximum_features"),4096);o.constant=requested.at("constant").boolean(true);
    o.normalize=requested.at("normalize").boolean();o.bandwidth=requested.at("bandwidth").number(1);
    return c::buildSpectralDictionary(data.x,data.dimension,o);
}
inline c::KoopmanSnapshotPairs pairs(const J& value) {
    c::KoopmanSnapshotPairs data;data.dimension=count(value.at("dimension"));data.x=numbers(value.at("x"));data.y=numbers(value.at("y"));
    if(!data.dimension||data.x.size()%data.dimension)throw std::invalid_argument("Invalid snapshot coordinates");
    data.samples=data.x.size()/data.dimension;data.weights=numbers(value.at("weights"));
    data.successorReplicates=numbers(value.at("successor_replicates"));data.replicateCount=count(value.at("replicate_count"));return data;
}
inline void measures(c::SpectralKROM& model,const J& settings) {
    for(const auto& v:settings.at("rigged_measures").array())
        c::riggedSpectralMeasure(model,numbers(v.at("coefficients")),v.at("name").string(),
            numbers(v.at("angles")),v.at("bandwidth").number(),count(v.at("order"),1));
}
} // namespace dynsyslab::workflow::spectral_settings
