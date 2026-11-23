#include <dynsyslab/math/Distributions.hpp>
#include <cmath>
#include <limits>
#include <stdexcept>
// Distribution differentiation and smooth multiplication: Hormander,
// The Analysis of Linear Partial Differential Operators I, chapter 2.
// Coarea: Federer, Geometric Measure Theory, 3.2.12. Fourier orbit pairing
// is elementary circle integration; no product of singular distributions.
namespace dynsyslab::math {
namespace {
void finite(std::complex<double> c){if(!std::isfinite(c.real())||!std::isfinite(c.imag()))throw std::invalid_argument("Finite distribution coefficients/test values required");}
void validate(const PointDistribution& d){
    if(!std::isfinite(d.location)||d.coefficients.size()>65)throw std::invalid_argument("Finite location and distribution order<=64 required");
    for(auto c:d.coefficients)finite(c);
}
void jet(const std::vector<std::complex<double>>& v,std::size_t n){if(v.size()<n)throw std::invalid_argument("Insufficient test/smooth derivatives");for(auto c:v)finite(c);}
}
std::complex<double> pairTestJet(const PointDistribution& d,const std::vector<std::complex<double>>& test){
    validate(d);jet(test,d.coefficients.size());std::complex<double> out{};
    for(std::size_t j=0;j<d.coefficients.size();++j)out+=(j%2?-1.:1.)*d.coefficients[j]*test[j];finite(out);return out;
}
PointDistribution distributionDerivative(const PointDistribution& d,unsigned order){
    validate(d);if(order>64||(!d.coefficients.empty()&&d.coefficients.size()+order>65))throw std::invalid_argument("Distribution derivative exceeds order 64");
    PointDistribution out{d.location,{}};if(d.coefficients.empty())return out;
    out.coefficients.resize(order);out.coefficients.insert(out.coefficients.end(),d.coefficients.begin(),d.coefficients.end());return out;
}
PointDistribution multiplySmooth(const PointDistribution& d,const std::vector<std::complex<double>>& smooth){
    validate(d);jet(smooth,d.coefficients.size());PointDistribution out{d.location,std::vector<std::complex<double>>(d.coefficients.size())};
    for(std::size_t j=0;j<d.coefficients.size();++j){double choose=1;
        for(std::size_t k=0;k<=j;++k){out.coefficients[j-k]+=(k%2?-1.:1.)*choose*d.coefficients[j]*smooth[k];
            if(k<j)choose*=double(j-k)/double(k+1);}}
    validate(out);return out;
}
std::complex<double> periodicOrbitPairing(const std::vector<std::complex<double>>& samples,int n,double period,bool normalized){
    if(n<-1000000||n>1000000||samples.empty()||std::size_t(std::abs(n))*2>=samples.size()||!std::isfinite(period)||period<=0)
        throw std::invalid_argument("Orbit pairing requires a positive period and unaliased uniformly phased samples");
    std::complex<double> out{};const double angle=2*std::acos(-1.)/double(samples.size());
    for(std::size_t j=0;j<samples.size();++j){finite(samples[j]);out+=std::exp(std::complex<double>(0,n*angle*j))*samples[j];}
    out*= (normalized?1.:period)/double(samples.size());finite(out);return out;
}
std::complex<double> levelSetDeltaPairing(const std::vector<std::complex<double>>& values,const std::vector<double>& weights,const std::vector<double>& norms){
    if(values.size()!=weights.size()||values.size()!=norms.size())throw std::invalid_argument("Coarea quadrature arrays must agree");
    std::complex<double> out{};
    for(std::size_t i=0;i<values.size();++i){finite(values[i]);if(!std::isfinite(weights[i])||weights[i]<0||!std::isfinite(norms[i])||norms[i]<=0)
        throw std::invalid_argument("Regular level-set quadrature requires nonnegative weights and positive gradient norms");out+=values[i]*weights[i]/norms[i];}
    finite(out);return out;
}
double spectralPushforwardDensity(const std::vector<SpectralDensityRoot>& roots){
    double out=0;
    for(const auto& r:roots){finite(r.observableCoefficient);
        if(!r.harmonic||!std::isfinite(r.frequencyDerivative)||r.frequencyDerivative==0||!std::isfinite(r.sourceDensity)||r.sourceDensity<0)
            throw std::invalid_argument("Regular nonzero-frequency roots required; atoms and critical roots need separate treatment");
        out+=r.sourceDensity*std::norm(r.observableCoefficient)/std::abs(double(r.harmonic)*r.frequencyDerivative);}
    if(!std::isfinite(out))throw std::overflow_error("Spectral density overflow");return out;
}
double gaussianDeltaDerivative(double x,double epsilon,unsigned order){
    if(!std::isfinite(x)||!std::isfinite(epsilon)||epsilon<=0||order>32)throw std::invalid_argument("Finite x, epsilon>0 and order<=32 required");
    const double z=x/epsilon;
    // For binary64 epsilon and order<=32, |z|>1000 underflows even after
    // epsilon^(-33) scaling. A smaller cutoff independent of epsilon is unsafe.
    if(!std::isfinite(z)||std::abs(z)>1000)return 0;
    double h0=1,h1=z;
    for(unsigned j=1;j<order;++j){const double next=z*h1-j*h0;h0=h1;h1=next;}
    const double hermite=order?h1:1;
    if(hermite==0)return 0;
    const double logMagnitude=std::log(std::abs(hermite))-z*z/2
        -.5*std::log(2*std::acos(-1.))-(order+1)*std::log(epsilon);
    if(logMagnitude>std::log(std::numeric_limits<double>::max()))
        throw std::overflow_error("Gaussian regularization overflow");
    return (order%2?-1.:1.)*std::copysign(std::exp(logMagnitude),hermite);
}
}
