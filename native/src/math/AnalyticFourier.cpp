#include <dynsyslab/math/AnalyticFourier.hpp>
#include <dynsyslab/core/SpecialFunctions.hpp>
#include <cmath>
#include <stdexcept>

// Jacobi Fourier series: DLMF 22.11.1-3, https://dlmf.nist.gov/22.11.
// Logarithmic series: log(1+rho*z), |rho|<1, plus its conjugate on |z|=1.
// Geometric majorants provide bounds for ALL omitted coefficients, not just
// a comparison with a larger finite truncation.
namespace dynsyslab::math {
namespace {
constexpr double pi=3.14159265358979323846;
void valid(int n,double w){if(n<0||n>1000000||!std::isfinite(w)||w<=0)
    throw std::invalid_argument("Need maximum harmonic 0..1000000 and positive finite frequency");}
void pair(FourierSeries& s,int n,std::complex<double> c){s.coefficients[n]=c;s.coefficients[-n]=std::conj(c);}
}
AnalyticFourierExpansion jacobiFourier(const std::string& kind,double m,int N,double w){
    valid(N,w);if((kind!="sn"&&kind!="cn"&&kind!="dn")||!std::isfinite(m)||m<0||m>=1)
        throw std::invalid_argument("Jacobi Fourier needs sn/cn/dn and 0 <= parameter < 1");
    const double K=core::special::ellipticK(m);
    AnalyticFourierExpansion out{{w,{}},0,kind=="dn"?2*K:4*K};
    if(m==0){if(kind=="dn")out.series.coefficients[0]=1;
        else if(N>=1)pair(out.series,1,kind=="sn"?std::complex<double>(0,-.5):.5);
        else out.uniformTailBound=1;return out;}
    // Complementary K by AGM(1,sqrt(m)); forming 1-m loses tiny m.
    // DLMF 19.8.5. This is the same complete integral, not an asymptotic fit.
    double a=1,b=std::sqrt(m);
    for(unsigned j=0;j<64;++j){const double next=(a+b)/2;b=std::sqrt(a)*std::sqrt(b);a=next;
        if(std::abs(a-b)<=2e-16*a)break;}
    const double Kprime=pi/(a+b);
    const double logq=-pi*Kprime/K;
    if(kind=="dn"){
        out.series.coefficients[0]=pi/(2*K);
        for(int n=1;n<=N;++n)pair(out.series,n,pi/K*std::exp(n*logq)/(1+std::exp(2*n*logq)));
        out.uniformTailBound=2*pi/K*std::exp((N+1)*logq)/(-std::expm1(logq));
    }else{
        const double scale=pi/(K*std::sqrt(m));
        for(int n=1;n<=N;n+=2){const double c=scale*std::exp(.5*n*logq)/
            (kind=="sn"?-std::expm1(n*logq):1+std::exp(n*logq));
            pair(out.series,n,kind=="sn"?std::complex<double>(0,-c):c);}
        const int next=(N%2)?N+2:N+1;
        out.uniformTailBound=2*scale*std::exp(.5*next*logq)/(-std::expm1(logq));
        if(kind=="sn")out.uniformTailBound/=(-std::expm1(next*logq));
    }
    if(!std::isfinite(out.uniformTailBound))throw std::overflow_error("Fourier tail bound overflow");return out;
}
AnalyticFourierExpansion logarithmicCosineFourier(double k,int N,double w){
    valid(N,w);if(!std::isfinite(k)||std::abs(k)>=1)throw std::invalid_argument("Need |k| < 1");
    const double s=std::sqrt((1-k)*(1+k)),rho=k/(1+s),a=std::abs(rho);
    AnalyticFourierExpansion out{{w,{{0,std::log1p((s-1)/2)}}},0,2*pi};
    double power=rho;for(int n=1;n<=N;++n){pair(out.series,n,(n%2?1.:-1.)*power/n);power*=rho;}
    if(a>0)out.uniformTailBound=2*std::pow(a,N+1)/((N+1)*(1-a));return out;
}
}
