#include <dynsyslab/core/PeriodicEstimation.hpp>
#include <dynsyslab/core/EnsembleBayes.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace dynsyslab::core {
PeriodicFilterResult filterPeriodicAtlas(const math::PeriodicAtlas& atlas,
    const std::vector<double>& times,const std::vector<double>& measurements,const PeriodicFilterOptions& o){
    const auto np=o.parameters.size(),nm=o.measurementChannels.size(),nc=atlas.knots()[0].channels.size();
    if(!np||!nm||o.outputChannels.empty()||times.empty()||measurements.size()!=times.size()*nm||
       o.measurementVariances.size()!=nm||o.phaseCount<3||o.tableSamples<3||
       np>1000000/o.phaseCount||np>10000000/o.tableSamples/nc)
        throw std::invalid_argument("Invalid periodic-filter dimensions or resource limit");
    for(auto c:o.measurementChannels)if(c>=nc)throw std::invalid_argument("Measurement channel out of range");
    for(auto c:o.outputChannels)if(c>=nc)throw std::invalid_argument("Output channel out of range");
    if(!std::isfinite(o.frozenFrequency)||o.frozenFrequency<0)throw std::invalid_argument("Frozen frequency must be nonnegative");
    for(std::size_t i=0;i<times.size();++i)if(!std::isfinite(times[i])||(i&&times[i]<=times[i-1]))throw std::invalid_argument("Filter times must increase strictly");
    for(std::size_t i=0;i<np;++i)if(i&&o.parameters[i]<=o.parameters[i-1])throw std::invalid_argument("Filter parameters must increase strictly");
    auto prior=o.parameterPrior;if(prior.empty())prior.assign(np,1.);
    if(prior.size()!=np)throw std::invalid_argument("Parameter prior size mismatch");
    const double twoPi=2*std::acos(-1.);std::vector<double> table(np*nc*o.tableSamples),omega;
    std::vector<std::vector<double>> particles;std::vector<double> particlePrior;
    for(std::size_t p=0;p<np;++p){const auto channels=atlas.interpolate(o.parameters[p]);
        omega.push_back(o.frozenFrequency>0?o.frozenFrequency:channels[0].angularFrequency);
        for(std::size_t c=0;c<nc;++c){const auto samples=math::uniformFourierSamples(channels[c],o.tableSamples);
            double scale=0,imaginary=0;for(auto z:samples){scale=std::max(scale,std::abs(z.real()));imaginary=std::max(imaginary,std::abs(z.imag()));}
            if(imaginary>1e-9*(1+scale))throw std::invalid_argument("Periodic filter needs real measurement/output channels");
            for(std::size_t j=0;j<samples.size();++j)table[(p*nc+c)*o.tableSamples+j]=samples[j].real();
        }
        for(std::size_t j=0;j<o.phaseCount;++j){particles.push_back({double(p),twoPi*j/o.phaseCount});particlePrior.push_back(prior[p]/o.phaseCount);}
    }
    EnsembleBayesFilter filter(std::move(particles),std::move(particlePrior));
    double elapsed=0;
    auto read=[&](const std::vector<double>& x,const std::vector<std::size_t>& channels){
        const auto p=std::size_t(x[0]);double phase=std::fmod(x[1]+omega[p]*elapsed,twoPi);if(phase<0)phase+=twoPi;
        const double cell=phase/twoPi*o.tableSamples;const auto i=std::min(o.tableSamples-1,std::size_t(cell)),j=(i+1)%o.tableSamples;
        const double a=cell-i;std::vector<double> out;
        for(auto c:channels){const auto base=(p*nc+c)*o.tableSamples;out.push_back((1-a)*table[base+i]+a*table[base+j]);}return out;};
    PeriodicFilterResult result;
    for(std::size_t k=0;k<times.size();++k){elapsed=times[k]-times.front();
        std::vector<double> y(measurements.begin()+k*nm,measurements.begin()+(k+1)*nm);
        result.logEvidence.push_back(filter.correct(y,o.measurementVariances,[&](const auto& x){return read(x,o.measurementChannels);}));
        const auto out=filter.mean([&](const auto& x){return read(x,o.outputChannels);});result.outputs.insert(result.outputs.end(),out.begin(),out.end());
        const auto mean=filter.mean([&](const auto& x){double p=o.parameters[std::size_t(x[0])];return std::vector<double>{p,p*p};});
        result.parameterMean.push_back(mean[0]);result.parameterVariance.push_back(std::max(0.,mean[1]-mean[0]*mean[0]));
        result.effectiveSampleSize.push_back(filter.effectiveSampleSize());
    }
    result.finalWeights=filter.weights();return result;
}
}
