#include "KoopmanSpectralInternal.hpp"
namespace dynsyslab::core {
namespace {
void validateReplicates(const KoopmanSnapshotPairs& data) {
    if(data.replicateCount<2||data.successorReplicates.size()!=data.replicateCount*data.samples*data.dimension)
        throw std::invalid_argument("Stochastic analysis needs at least two independent successors per source snapshot");
}
}
KoopmanModel fitStochasticKoopmanSpectral(const KoopmanSnapshotPairs& data,const SpectralDictionary& dictionary,
                                          const KoopmanSpectralOptions& options,std::optional<FunctionDomain> domain) {
    spectral_detail::check(data,dictionary.dimension);validateReplicates(data);
    KoopmanSnapshotPairs pooled;pooled.samples=data.samples*data.replicateCount;pooled.dimension=data.dimension;
    const auto w=spectral_detail::weights(data.samples,data.weights,options.weightPolicy);
    pooled.y=data.successorReplicates;
    for(std::size_t j=0;j<data.replicateCount;++j) {
        pooled.x.insert(pooled.x.end(),data.x.begin(),data.x.end());
        for(double v:w)pooled.weights.push_back(v/data.replicateCount);
    }
    auto settings=options;settings.method="resdmd";settings.weightPolicy="uniform";
    auto model=fitKoopmanSpectral(pooled,dictionary,settings,std::move(domain));
    auto compression=*model.compression();compression.method="stochastic_resdmd";
    compression.assumptions="Replicated independent successor realizations at identical source points asserted by caller. Regression approximates conditional expectation.";
    model.setCompression(std::move(compression));
    auto weightedData=data;weightedData.weights=w;
    auto diagnostics=model.spectralDiagnostics();diagnostics.stochastic=stochasticKoopmanResiduals(model,weightedData);
    model.setSpectralDiagnostics(std::move(diagnostics));
    return model;
}
StochasticSpectralResiduals stochasticKoopmanResiduals(const KoopmanModel& model,const KoopmanSnapshotPairs& data) {
    if(!model.compression())throw std::invalid_argument("Stochastic residuals need the krom's compressed operator");
    const auto& dictionary=model.compression()->dictionary;
    spectral_detail::check(data,dictionary.dimension);validateReplicates(data);
    // The krom's eigenfunctions on its frozen dictionary, at the sources and at every successor.
    const auto& coefficients=frozenDictionaryExpansion(model).coefficients();
    const auto multipliers=model.spectrum().discreteMultipliers();
    const auto r=multipliers.size(),n=dictionary.size();
    const auto evaluate=[&](const std::vector<double>& points) {
        const auto phi=dictionary.evaluate(points);const auto count=points.size()/dictionary.dimension;
        std::vector<std::complex<double>> out(count*r);
        for(std::size_t i=0;i<count;++i)for(std::size_t j=0;j<r;++j)for(std::size_t k=0;k<n;++k)out[i*r+j]+=phi[i*n+k]*coefficients[j*n+k];
        return out;
    };
    const auto x=evaluate(data.x),y=evaluate(data.successorReplicates);
    const auto w=spectral_detail::weights(data.samples,data.weights);
    const auto m=data.samples,R=data.replicateCount;
    StochasticSpectralResiduals out;out.replicates=R;
    out.interpretation="Sample-mean residual and unbiased innovation standard deviation are distinct. Total is RMS pathwise residual; total^2=sampleMean^2+(R-1)/R*innovation^2. Mean residual retains finite-replicate sampling error.";
    for(std::size_t j=0;j<r;++j) {
        double meanError=0,variance=0,denom=0;
        for(std::size_t i=0;i<m;++i) {
            std::complex<double> mean{};
            for(std::size_t k=0;k<R;++k)mean+=y[(k*m+i)*r+j]/double(R);
            meanError+=w[i]*std::norm(mean-multipliers[j]*x[i*r+j]);denom+=w[i]*std::norm(x[i*r+j]);
            for(std::size_t k=0;k<R;++k)variance+=w[i]*std::norm(y[(k*m+i)*r+j]-mean)/(R-1);
        }
        if(denom<=1e-28)throw std::invalid_argument("Stochastic residual denominator vanishes");
        out.sampleMean.push_back(std::sqrt(meanError/denom));out.innovation.push_back(std::sqrt(variance/denom));
        out.total.push_back(std::sqrt((meanError+double(R-1)/R*variance)/denom));
    }
    return out;
}
}
