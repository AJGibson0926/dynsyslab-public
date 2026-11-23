#include <dynsyslab/core/CoefficientShift.hpp>
#include <dynsyslab/core/DiscreteFourier.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>
// Fourier coefficient covariance under time translation; Koopman projection:
// Mezic, Nonlinear Dynamics 41 (2005), https://doi.org/10.1007/s11071-005-2824-x.
namespace dynsyslab::core {
namespace {
std::vector<std::vector<std::complex<double>>> coefficients(const PeriodicSamples& x,int maximum) {
    const auto count=x.times.size();
    std::vector<std::vector<std::complex<double>>> result(x.channels);
    for(std::size_t channel=0;channel<x.channels;++channel) {
        std::vector<std::complex<double>> trace(count);
        for(std::size_t i=0;i<count;++i)trace[i]=x.observables[i*x.channels+channel];
        const auto spectrum=discreteFourierTransform(trace);
        for(int n=-maximum;n<=maximum;++n)
            result[channel].push_back(spectrum[n<0?count+n:std::size_t(n)]/double(count));
    }return result;
}
}
CoefficientShiftCheck checkPeriodicCoefficientShift(const AnalyticSolutionFamily& family,
    const AnalyticState& initial,const std::vector<Observable>& obs,std::size_t count,
    int maximum,const std::vector<double>& shifts,double t0,double explicitPeriod) {
    if(initial.empty()||maximum<0||count<3||2*std::size_t(maximum)>=count||
       shifts.empty()||shifts.size()>256||count>1000000||obs.empty()||
       count>100000000/(shifts.size()+1)/initial.size())
        throw std::invalid_argument("Invalid or oversized coefficient-shift request");
    for(double x:initial)if(!std::isfinite(x))throw std::invalid_argument("Finite IC required");
    for(double s:shifts)if(!std::isfinite(s))throw std::invalid_argument("Finite shifts required");
    const auto base=family.instantiate(initial,t0);
    if(!base)throw std::invalid_argument("Inapplicable analytic family: "+base.reason);
    const auto sampled=samplePeriodicTrajectory(base,initial.size(),obs,count,t0,explicitPeriod);
    if(sampled.channels>100000000/(shifts.size()+1)/count)
        throw std::invalid_argument("Coefficient-shift measurement count exceeds limit");
    const auto reference=coefficients(sampled,maximum);
    CoefficientShiftCheck out;out.period=sampled.period;out.angularFrequency=sampled.angularFrequency;
    out.closureError=sampled.closureError;out.channels=sampled.channels;out.samples=count;out.maximumHarmonic=maximum;
    for(double shift:shifts) {
        const auto state=base.evaluate(t0+shift);
        if(state.size()!=initial.size())throw std::domain_error("Shifted state dimension changed");
        for(double x:state)if(!std::isfinite(x))throw std::domain_error("Nonfinite shifted state");
        const auto restarted=family.instantiate(state,t0);
        if(!restarted)throw std::domain_error("Shifted IC inapplicable: "+restarted.reason);
        const auto target=samplePeriodicTrajectory(restarted,initial.size(),obs,count,t0,explicitPeriod);
        const auto actual=coefficients(target,maximum);
        CoefficientShiftSample row;row.shift=shift;row.closureError=target.closureError;
        row.relativeFrequencyError=std::abs(target.angularFrequency-sampled.angularFrequency)/sampled.angularFrequency;
        for(std::size_t j=0;j<state.size();++j)
            row.initialStateError=std::max(row.initialStateError,std::abs(state[j]-target.states[j])/(1+std::abs(state[j])));
        for(std::size_t channel=0;channel<out.channels;++channel) {
            double error=0,norm=0,maximumError=0;
            for(int n=-maximum;n<=maximum;++n) {
                const auto z=reference[channel][n+maximum];
                const double phase=n*sampled.angularFrequency*shift;
                if(!std::isfinite(phase))throw std::overflow_error("Coefficient-shift phase overflow");
                const auto expected=z*std::exp(std::complex<double>(0,std::remainder(phase,2*std::acos(-1.))));
                const double e=std::abs(actual[channel][n+maximum]-expected);
                error=std::hypot(error,e);norm=std::hypot(norm,std::abs(z));maximumError=std::max(maximumError,e);
            }
            row.maximumCoefficientError.push_back(maximumError);row.coefficientL2Error.push_back(error);row.referenceL2Norm.push_back(norm);
        }out.shifts.push_back(std::move(row));
    }return out;
}
}
