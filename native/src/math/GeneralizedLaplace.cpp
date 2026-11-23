#include <dynsyslab/math/GeneralizedLaplace.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
namespace dynsyslab::math {
namespace {
using C=std::complex<double>;
const double pi=std::acos(-1.);
void validate(const SampledSignal& g){
    if(g.times.size()<4||g.times.size()!=g.values.size())throw std::invalid_argument("Signal needs at least four matching samples");
    for(std::size_t k=1;k<g.times.size();++k)if(!std::isfinite(g.times[k])||g.times[k]<=g.times[k-1])throw std::invalid_argument("Times must increase strictly");
    if(!std::isfinite(g.times[0]))throw std::invalid_argument("Times must be finite");
    for(const auto& v:g.values)if(!std::isfinite(v.real())||!std::isfinite(v.imag()))throw std::invalid_argument("Signal values must be finite");
}
double bump(double x){return x<=0||x>=1?0:std::exp(-1/(x*(1-x)));}
// Trapezoid weights on the nonuniform grid, times the bump over [a,b].
std::vector<double> windowWeights(const std::vector<double>& t,double a,double b){
    std::vector<double> w(t.size());double total=0;
    for(std::size_t k=0;k<t.size();++k){
        const double left=k?t[k]-t[k-1]:0,right=k+1<t.size()?t[k+1]-t[k]:0;
        w[k]=.5*(left+right)*bump((t[k]-a)/(b-a));total+=w[k];
    }
    if(!(total>0))throw std::invalid_argument("Averaging window contains no samples");
    for(auto& x:w)x/=total;
    return w;
}
C weightedExponentialMean(const std::vector<double>& t,const ComplexVector& y,const std::vector<double>& w,C lambda,double origin){
    C sum{};for(std::size_t k=0;k<t.size();++k)if(w[k]>0)sum+=w[k]*y[k]*std::exp(-lambda*(t[k]-origin));
    return sum;
}
double rms(const std::vector<double>& t,const ComplexVector& y){
    double sum=0;for(std::size_t k=0;k+1<t.size();++k)sum+=.5*(t[k+1]-t[k])*(std::norm(y[k])+std::norm(y[k+1]));
    return std::sqrt(sum/(t.back()-t.front()));
}
}
GeneralizedLaplaceResult generalizedLaplaceAnalysis(const SampledSignal& g,ComplexVector eigenvalues,const GeneralizedLaplaceOptions& o){
    validate(g);if(eigenvalues.empty())throw std::invalid_argument("GLA needs eigenvalues");
    if(!(o.amplificationLimit>1)||!std::isfinite(o.amplificationLimit))throw std::invalid_argument("Amplification limit must exceed one");
    for(const auto& z:eigenvalues)if(!std::isfinite(z.real())||!std::isfinite(z.imag()))throw std::invalid_argument("Eigenvalues must be finite");
    // Slow to fast; ties by |Im| for a reproducible order.
    std::stable_sort(eigenvalues.begin(),eigenvalues.end(),[](C a,C b){
        return a.real()!=b.real()?a.real()>b.real():std::abs(a.imag())<std::abs(b.imag());});
    for(std::size_t j=1;j<eigenvalues.size();++j)if(std::abs(eigenvalues[j]-eigenvalues[j-1])<1e-14*(1+std::abs(eigenvalues[j])))throw std::invalid_argument("GLA eigenvalues must be distinct");
    const double t0=g.times.front(),span=g.times.back()-t0;
    const std::size_t M=eigenvalues.size(),N=g.times.size();
    std::vector<std::vector<double>> weights(M);GeneralizedLaplaceResult out;out.modes.resize(M);
    for(std::size_t j=0;j<M;++j){
        const double rate=eigenvalues[j].real();
        double horizon=span;
        if(rate<0)horizon=std::min(span,std::log(o.amplificationLimit)/(-rate));
        // Keep several samples inside even a fast mode's short window.
        horizon=std::max(horizon,span*std::min(1.,16./double(N)));
        weights[j]=windowWeights(g.times,t0,t0+horizon);
        out.modes[j].eigenvalue=eigenvalues[j];out.modes[j].horizon=horizon;
    }
    ComplexVector residual=g.values;
    const auto subtract=[&](std::size_t j,C amount){
        if(amount==C{})return;
        for(std::size_t k=0;k<N;++k)residual[k]-=amount*std::exp(eigenvalues[j]*(g.times[k]-t0));
    };
    // Classical GLA: one ordered pass, each mode averaged after subtracting the slower ones.
    for(std::size_t j=0;j<M;++j){
        const C c=weightedExponentialMean(g.times,residual,weights[j],eigenvalues[j],t0);
        out.modes[j].coefficient=out.modes[j].sequentialCoefficient=c;subtract(j,c);
    }
    if(o.refine){
        // The finite-window averages obey G c=b with G_ji=mean_j[e^{(lambda_i-lambda_j)(t-t0)}].
        // Solving this directly removes cross-leakage; plain re-sweeps of it can diverge.
        ComplexMatrix G(M,M);ComplexVector b(M);
        for(std::size_t j=0;j<M;++j){
            b[j]=weightedExponentialMean(g.times,g.values,weights[j],eigenvalues[j],t0);
            for(std::size_t i=0;i<M;++i){C sum{};
                for(std::size_t k=0;k<N;++k)if(weights[j][k]>0)sum+=weights[j][k]*std::exp((eigenvalues[i]-eigenvalues[j])*(g.times[k]-t0));
                G(j,i)=sum;}
            // Row equilibration: e^{-lambda_j t} makes the rows of fast modes very large.
            double scale=0;for(std::size_t i=0;i<M;++i)scale=std::max(scale,std::abs(G(j,i)));
            if(scale>0){for(std::size_t i=0;i<M;++i)G(j,i)/=scale;b[j]/=scale;}
        }
        const auto sv=complexSingularValues(G);
        out.refinementCondition=sv.back()>0?sv.front()/sv.back():std::numeric_limits<double>::infinity();
        if(out.refinementCondition<o.maximumCondition){
            const auto c=ComplexLu(G).solve(b);out.refined=true;residual=g.values;
            for(std::size_t j=0;j<M;++j){out.modes[j].coefficient=c[j];subtract(j,c[j]);}
        }
    }
    out.residual=std::move(residual);out.residualRms=rms(g.times,out.residual);out.signalRms=rms(g.times,g.values);
    return out;
}
std::map<int,C> estimateHarmonics(const SampledSignal& g,double w,int maximumHarmonic,double fraction){
    validate(g);
    if(!(w>0)||!std::isfinite(w))throw std::invalid_argument("Frequency must be positive");
    if(maximumHarmonic<0||maximumHarmonic>4096)throw std::invalid_argument("Harmonic bound must be 0..4096");
    if(!(fraction>=0&&fraction<1))throw std::invalid_argument("Transient fraction must lie in [0,1)");
    const double T=2*pi/w,end=g.times.back(),raw=g.times.front()+fraction*(end-g.times.front());
    const double periods=std::floor((end-raw)/T);
    if(periods<2)throw std::invalid_argument("Tail must hold at least two whole periods");
    const auto weight=windowWeights(g.times,end-periods*T,end);
    std::map<int,C> out;
    for(int n=-maximumHarmonic;n<=maximumHarmonic;++n)out[n]=weightedExponentialMean(g.times,g.values,weight,C(0,n*w),0);
    return out;
}
StationaryDecomposition decomposeStationary(const SampledSignal& g,double w,int maximumHarmonic,double fraction,double floor){
    StationaryDecomposition out;out.harmonics=estimateHarmonics(g,w,maximumHarmonic,fraction);
    if(!(floor>=0)||!std::isfinite(floor))throw std::invalid_argument("Noise floor must be nonnegative");
    const std::size_t N=g.times.size();out.stationary.resize(N);out.transient.resize(N);
    double scale=0;
    for(std::size_t k=0;k<N;++k){C sum{};for(const auto& [n,c]:out.harmonics)sum+=c*std::exp(C(0,n*w*g.times[k]));
        out.stationary[k]=sum;out.transient[k]=g.values[k]-sum;scale=std::max(scale,std::abs(out.transient[k]));}
    // Stroboscopic samples t0+kT by linear interpolation of the complex transient.
    const double T=2*pi/w,t0=g.times.front();std::vector<double> xs,ys;
    // The tail fixes the harmonics, so |transient| there is their error plateau, not decay.
    const double tailStart=t0+fraction*(g.times.back()-t0);double plateau=0;
    for(std::size_t k=0;k<N;++k)if(g.times[k]>=tailStart)plateau=std::max(plateau,std::abs(out.transient[k]));
    const double threshold=std::max(floor*scale,10*plateau);
    std::size_t cursor=0;
    for(std::size_t k=0;;++k){const double t=t0+double(k)*T;if(t>g.times.back())break;
        while(cursor+2<N&&g.times[cursor+1]<t)++cursor;
        const double h=g.times[cursor+1]-g.times[cursor],x=(t-g.times[cursor])/h;
        const double a=std::abs((1-x)*out.transient[cursor]+x*out.transient[cursor+1]);
        if(t<tailStart&&a>threshold&&a>0){xs.push_back(t);ys.push_back(std::log(a));}
    }
    out.strobeCount=xs.size();
    if(xs.size()>=3){
        const double mx=std::accumulate(xs.begin(),xs.end(),0.)/double(xs.size()),my=std::accumulate(ys.begin(),ys.end(),0.)/double(ys.size());
        double sxx=0,sxy=0;for(std::size_t k=0;k<xs.size();++k){sxx+=(xs[k]-mx)*(xs[k]-mx);sxy+=(xs[k]-mx)*(ys[k]-my);}
        out.decayRate=sxy/sxx;double r=0;
        for(std::size_t k=0;k<xs.size();++k){const double e=ys[k]-my-out.decayRate*(xs[k]-mx);r+=e*e;}
        out.fitResidual=std::sqrt(r/double(xs.size()));out.decayRateAvailable=true;
    }
    return out;
}
PoleOrderProbe probePoleOrder(const SampledSignal& g,C lambda){
    validate(g);const double a=g.times.front(),b=g.times.back(),m=.5*(a+b);
    const auto first=weightedExponentialMean(g.times,g.values,windowWeights(g.times,a,m),lambda,0);
    const auto second=weightedExponentialMean(g.times,g.values,windowWeights(g.times,m,b),lambda,0);
    PoleOrderProbe out;out.ratio=std::abs(first)>0?std::abs(second)/std::abs(first):0;
    out.order=out.ratio>2?2:1;
    return out;
}
}
