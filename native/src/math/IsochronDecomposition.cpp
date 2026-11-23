#include <dynsyslab/math/IsochronDecomposition.hpp>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
namespace dynsyslab::math {
namespace {
using C=std::complex<double>;
const double pi=std::acos(-1.);
// Harmonic with the largest |g_m|, m != 0: the best-conditioned phase reference.
int referenceHarmonic(const FourierSeries& s){
    int best=0;double top=0;
    for(const auto& [m,g]:s.coefficients){if(m==0)continue;const double a=std::abs(g);
        if(a>top||(a==top&&m>0&&best<0)){top=a;best=m;}}
    return best;
}
// The series may be invariant under theta -> theta+2 pi/d when all harmonics share a factor d.
int harmonicGcd(const FourierSeries& s){
    int d=0;for(const auto& [m,g]:s.coefficients){if(m==0||g==C{})continue;int a=std::abs(m),b=d;while(b){const int r=a%b;a=b;b=r;}d=a;}
    return d;
}
}
IsochronSplit isochronSplit(const SampledSignal& g,const FourierSeries& cycle,const IsochronSplitOptions& o){
    if(g.times.size()<8||g.times.size()!=g.values.size())throw std::invalid_argument("Isochron split needs at least eight matching samples");
    const double W=cycle.angularFrequency;
    if(!(W>0)||!std::isfinite(W))throw std::invalid_argument("Cycle frequency must be positive");
    if(cycle.coefficients.empty())throw std::invalid_argument("Cycle needs Fourier coefficients");
    if(!(o.noiseFloor>=0)||!(o.phaseTolerance>0)||!(o.transientFraction>=0&&o.transientFraction<1))throw std::invalid_argument("Invalid isochron options");
    if(o.asymptoticPhase&&!std::isfinite(*o.asymptoticPhase))throw std::invalid_argument("Asymptotic phase must be finite");
    const int m=referenceHarmonic(cycle),symmetry=harmonicGcd(cycle);
    if(m==0||symmetry==0)throw std::invalid_argument("Cycle needs a nonzero harmonic to define a phase");
    const double t0=g.times.front();IsochronSplit out;
    // Tail estimate: c_m e^{i m W (t-t0)} against g_m e^{i m theta0}.
    SampledSignal shifted=g;for(auto& t:shifted.times)t-=t0;shifted.derivatives.clear();
    try{
        const auto tail=estimateHarmonics(shifted,W,std::abs(m),o.transientFraction);
        const C ratio=tail.at(m)/cycle.coefficients.at(m);
        if(std::abs(ratio)>1e-6)out.tailPhase=std::arg(ratio)/m;
    }catch(const std::invalid_argument&){out.warnings.push_back("Tail holds fewer than two periods: no independent phase check");}
    if(o.asymptoticPhase){out.asymptoticPhase=*o.asymptoticPhase;out.phaseSource="phase_map";}
    else if(out.tailPhase){out.asymptoticPhase=*out.tailPhase;out.phaseSource="tail_estimate";
        out.warnings.push_back("No phase map supplied: the phase is fitted from the same record, so the split is not an independent test");}
    else throw std::invalid_argument("No phase map and no usable tail to estimate the asymptotic phase");
    if(symmetry>1)out.warnings.push_back("Cycle harmonics share a factor: this observable fixes the phase only modulo 2*pi/"+std::to_string(symmetry));
    if(out.tailPhase){const double width=2*pi/symmetry;
        out.phaseMismatch=std::abs(std::remainder(out.asymptoticPhase-*out.tailPhase,width));
        if(*out.phaseMismatch>o.phaseTolerance)out.warnings.push_back("Supplied phase disagrees with the steady tail: wrong isochron, wrong cycle, or not converged");}
    const auto N=g.times.size();out.stationary.resize(N);out.transient.resize(N);double scale=0;
    for(std::size_t k=0;k<N;++k){out.stationary[k]=evaluate(cycle,out.asymptoticPhase+W*(g.times[k]-t0));
        out.transient[k]=g.values[k]-out.stationary[k];scale=std::max(scale,std::abs(out.transient[k]));}
    // Plateau: the final tenth holds the coefficient and phase error, not the decay.
    double plateau=0;for(std::size_t k=N-N/10;k<N;++k)plateau=std::max(plateau,std::abs(out.transient[k]));
    const double threshold=std::max(o.noiseFloor*scale,10*plateau),T=2*pi/W;std::vector<double> xs,ys;std::size_t cursor=0;
    for(std::size_t k=0;;++k){const double t=t0+double(k)*T;if(t>g.times.back())break;
        while(cursor+2<N&&g.times[cursor+1]<t)++cursor;
        const double x=(t-g.times[cursor])/(g.times[cursor+1]-g.times[cursor]);
        const double a=std::abs((1-x)*out.transient[cursor]+x*out.transient[cursor+1]);
        if(a>threshold&&a>0){xs.push_back(t);ys.push_back(std::log(a));}}
    out.strobeCount=xs.size();
    if(xs.size()>=3){
        const double mx=std::accumulate(xs.begin(),xs.end(),0.)/double(xs.size()),my=std::accumulate(ys.begin(),ys.end(),0.)/double(ys.size());
        double sxx=0,sxy=0;for(std::size_t k=0;k<xs.size();++k){sxx+=(xs[k]-mx)*(xs[k]-mx);sxy+=(xs[k]-mx)*(ys[k]-my);}
        const double rate=sxy/sxx;double r=0;for(std::size_t k=0;k<xs.size();++k){const double e=ys[k]-my-rate*(xs[k]-mx);r+=e*e;}
        out.decayRate=rate;out.fitResidual=std::sqrt(r/double(xs.size()));
        if(rate>=0)out.warnings.push_back("Transient does not decay: the orbit is not converging to this cycle at this phase");
    }else out.warnings.push_back("Too few stroboscopic samples above the error plateau for a decay rate");
    return out;
}
IsochronSplit isochronSplitFromMap(const SampledSignal& g,const FourierSeries& cycle,const PhaseMap& map,const std::vector<double>& x0,IsochronSplitOptions o){
    if(!map)throw std::invalid_argument("Phase map required");
    if(x0.empty()||!std::all_of(x0.begin(),x0.end(),[](double v){return std::isfinite(v);}))throw std::invalid_argument("Finite initial state required");
    o.asymptoticPhase=map(x0);
    return isochronSplit(g,cycle,o);
}
}
