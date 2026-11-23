#include <dynsyslab/core/PeriodicAnalysis.hpp>
#include <cmath>
#include <stdexcept>
namespace dynsyslab::core {
PeriodicSamples samplePeriodicTrajectory(const AnalyticTrajectory& trajectory,std::size_t d,
    const std::vector<Observable>& observables,std::size_t count,double t0,double explicitPeriod){
    if(!trajectory||!d||count<3||count>1000000||!std::isfinite(t0)||!std::isfinite(explicitPeriod)||explicitPeriod<0)
        throw std::invalid_argument("Invalid periodic trajectory sampling request");
    PeriodicSamples result;result.period=explicitPeriod>0?explicitPeriod:trajectory.period;
    if(!std::isfinite(result.period)||result.period<=0)throw std::invalid_argument("A known positive model-time period is required");
    result.angularFrequency=2*std::acos(-1.)/result.period;result.stateDimension=d;
    for(const auto& o:observables){if(o.stateDimension()!=d)throw std::invalid_argument("Periodic observable dimension mismatch");result.channels+=o.outputDimension();}
    if(result.channels>10000000/count||d>10000000/count)throw std::invalid_argument("Periodic sample output exceeds ten million entries");
    auto state=[&](double t){auto x=trajectory.evaluate(t);
        if(x.size()!=d)throw std::runtime_error("Periodic evaluator changed state dimension");
        for(double v:x)if(!std::isfinite(v))throw std::domain_error("Nonfinite periodic state");return x;};
    const auto first=state(t0),last=state(t0+result.period);
    for(std::size_t j=0;j<d;++j)result.closureError=std::max(result.closureError,std::abs(first[j]-last[j])/(1+std::abs(first[j])));
    for(std::size_t i=0;i<count;++i){
        const double t=t0+result.period*double(i)/count;auto x=state(t);result.times.push_back(t);
        result.states.insert(result.states.end(),x.begin(),x.end());
        for(const auto& o:observables){const auto values=o.evaluate(x.data());
            for(auto v:values){if(!std::isfinite(v.real())||!std::isfinite(v.imag()))throw std::domain_error("Nonfinite periodic observable");result.observables.push_back(v);}
        }
    }return result;
}
}
