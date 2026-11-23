#include <dynsyslab/math/PeriodicClock.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
// Change of variable in a Fourier coefficient, plus the Taylor expansion of
// exp(-i*n*beta). Standard Fourier algebra: Katznelson, Harmonic Analysis, ch. I.
// Sundman clocks in bubble solutions: Kudryashov & Sinelshchikov (2014),
// doi:10.1088/1751-8113/47/40/405202, equations (7),(25).
namespace dynsyslab::math {
namespace {
using C=std::complex<double>;
double l1(const FourierSeries& s){double a=0;for(const auto& [n,c]:s.coefficients)a+=std::abs(c);return a;}
struct Clock {double mean, nonconstantNorm, derivativeNorm;FourierSeries beta;};
Clock normalizedClock(const FourierSeries& w){
    meanSquare(w);const auto z=w.coefficients.find(0);
    if(z==w.coefficients.end()||z->second.imag()!=0||z->second.real()<=0)
        throw std::invalid_argument("Clock requires positive real mean");
    Clock c{z->second.real(),0,0,{w.angularFrequency,{}}};
    // Exact Hermitian symmetry makes the polynomial and beta real. Do not
    // silently discard imaginary input or repair an inconsistent supplied clock.
    for(const auto& [n,a]:w.coefficients)if(n){
        const auto opposite=w.coefficients.find(-n);
        if(opposite==w.coefficients.end()||opposite->second!=std::conj(a))
            throw std::invalid_argument("Clock Fourier coefficients must be Hermitian");
        c.nonconstantNorm+=std::abs(a);c.derivativeNorm+=std::abs(double(n)*a);
        c.beta.coefficients[n]=a/C(0,n*c.mean);
    }
    C zero{};for(const auto& [n,a]:c.beta.coefficients)zero-=a;
    // Summing conjugate pairs separately avoids an imaginary roundoff offset.
    c.beta.coefficients[0]=C(zero.real(),0);
    return c;
}
void compatible(const FourierSeries& g,const FourierSeries& w){
    meanSquare(g);if(g.angularFrequency!=w.angularFrequency)
        throw std::invalid_argument("Observable and clock auxiliary frequencies differ");
}
struct Sum {C value{},error{};void add(C x){const auto y=x-error,t=value+y;error=(t-value)-y;value=t;}};
}
PeriodicClockFourier reparameterizePeriodicFourier(const std::vector<FourierSeries>& channels,
    const FourierSeries& w,std::size_t samples,int maximum){
    if(channels.empty()||samples<3||samples>1000000||maximum<0||maximum>1000000||
       std::size_t(maximum)*2>=samples||channels.size()>100000000/samples/(2*std::size_t(maximum)+1))
        throw std::invalid_argument("Invalid periodic clock dimensions or work limit");
    const auto clock=normalizedClock(w);const auto rates=uniformFourierSamples(w,samples);
    const auto beta=uniformFourierSamples(clock.beta,samples);const double pi=std::acos(-1.);
    double minimum=std::numeric_limits<double>::infinity();for(auto r:rates)minimum=std::min(minimum,r.real());
    const double lower=std::max(clock.mean-clock.nonconstantNorm,minimum-pi/samples*clock.derivativeNorm);
    if(!(lower>0))throw std::domain_error("Supplied clock polynomial positivity not established; increase grid resolution or repair clock");
    PeriodicClockFourier out;out.meanClockRate=clock.mean;out.phaseCorrection=clock.beta;
    out.physicalPeriod=2*pi*clock.mean/w.angularFrequency;
    out.minimumSampledClockRate=minimum;out.clockPolynomialLowerBound=lower;
    for(const auto& g:channels){compatible(g,w);const auto values=uniformFourierSamples(g,samples);
        std::vector<Sum> coefficients(2*std::size_t(maximum)+1);
        for(std::size_t j=0;j<samples;++j){
            const double theta=2*pi*j/samples+beta[j].real();const auto phase=std::exp(C(0,-theta));
            const auto weighted=values[j]*(rates[j].real()/clock.mean/double(samples));
            coefficients[maximum].add(weighted);C power=1;
            for(int n=1;n<=maximum;++n){power*=phase;
                coefficients[maximum+n].add(weighted*power);
                coefficients[maximum-n].add(weighted*std::conj(power));}
        }
        FourierSeries s{w.angularFrequency/clock.mean,{}};
        for(int n=-maximum;n<=maximum;++n)s.coefficients[n]=coefficients[maximum+n].value;
        meanSquare(s);out.channels.push_back(std::move(s));
    }return out;
}
ClockCoefficientExpansion periodicClockCoefficientExpansion(const FourierSeries& g,
    const FourierSeries& w,int n,unsigned order){
    if(std::abs(double(n))>1000000||order>128)throw std::invalid_argument("Invalid clock expansion harmonic/order");
    const auto clock=normalizedClock(w);compatible(g,w);
    auto normalized=w;for(auto& [j,c]:normalized.coefficients)c/=clock.mean;
    auto term=multiply(g,normalized);Sum sum;
    for(unsigned k=0;k<=order;++k){const auto it=term.coefficients.find(n);if(it!=term.coefficients.end())sum.add(it->second);
        if(k<order&&n){term=multiply(term,clock.beta);for(auto& [j,c]:term.coefficients)c*=C(0,-double(n)/(k+1));}
        if(n==0)break;
    }
    const double scale=l1(g)*l1(normalized),b=std::abs(double(n))*l1(clock.beta);
    double bound=0;
    if(scale>0&&b>0){const double logBound=std::log(scale)+(order+1)*std::log(b)-std::lgamma(order+2.);
        if(logBound>std::log(std::numeric_limits<double>::max()))throw std::overflow_error("Clock Taylor bound overflows; use quadrature");
        bound=std::exp(logBound);}
    if(!std::isfinite(sum.value.real())||!std::isfinite(sum.value.imag()))throw std::overflow_error("Clock coefficient expansion overflow");
    return {sum.value,bound};
}
}
