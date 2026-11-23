#include <dynsyslab/math/HarmonicAnalysis.hpp>
#include <dynsyslab/core/DiscreteFourier.hpp>
#include <cmath>
#include <limits>
#include <stdexcept>
// Fourier series calculus and Parseval; Katznelson, An Introduction to
// Harmonic Analysis, chapter I. Uses the suite's existing FFT implementation.
namespace dynsyslab::math {
namespace {
void validate(const FourierSeries& s){
    if(!std::isfinite(s.angularFrequency)||s.angularFrequency<=0)throw std::invalid_argument("Positive finite angular frequency required");
    for(const auto& [n,c]:s.coefficients)if(n<-1000000||n>1000000||!std::isfinite(c.real())||!std::isfinite(c.imag()))
        throw std::invalid_argument("Invalid Fourier coefficient or harmonic index");
}
void finite(std::complex<double> c){if(!std::isfinite(c.real())||!std::isfinite(c.imag()))throw std::overflow_error("Fourier arithmetic overflow");}
}
std::complex<double> evaluate(const FourierSeries& s,double phase){
    validate(s);if(!std::isfinite(phase))throw std::invalid_argument("Finite phase required");
    phase=std::remainder(phase,2*std::acos(-1.));std::complex<double> out{};
    for(const auto& [n,c]:s.coefficients)out+=c*std::exp(std::complex<double>(0,n*phase));finite(out);return out;
}
std::vector<std::complex<double>> evolveFourierSeries(const FourierSeries& s,
    const std::vector<double>& times,double phase,double t0){
    validate(s);if(!std::isfinite(phase)||!std::isfinite(t0)||times.size()>1000000||
        (!s.coefficients.empty()&&times.size()>100000000/s.coefficients.size()))
        throw std::invalid_argument("Invalid or oversized Fourier evolution request");
    std::vector<std::complex<double>> out;out.reserve(times.size());
    for(double t:times){if(!std::isfinite(t))throw std::invalid_argument("Finite times required");
        out.push_back(evaluate(s,phase+s.angularFrequency*(t-t0)));}
    return out;
}
FourierSeries timeDerivative(const FourierSeries& s,unsigned order){
    validate(s);if(order>64)throw std::invalid_argument("Derivative order limited to 64");auto out=s;
    for(auto& [n,c]:out.coefficients){
        if(order&&n==0){c=0;continue;}
        const std::complex<double> factor(0,n*s.angularFrequency);
        for(unsigned j=0;j<order;++j)c*=factor;
        finite(c);
    }return out;
}
FourierSeries periodicPrimitive(const FourierSeries& s){
    validate(s);auto out=s;
    if(out.coefficients.count(0)&&out.coefficients.at(0)!=std::complex<double>{})
        throw std::invalid_argument("Nonzero mean has no periodic primitive; the secular term must be represented separately");
    for(auto& [n,c]:out.coefficients)if(n){
        const double frequency=n*s.angularFrequency;
        if(!std::isfinite(frequency))throw std::overflow_error("Harmonic frequency overflow");
        c={c.imag()/frequency,-c.real()/frequency};finite(c);
    }return out;
}
FourierSeries multiply(const FourierSeries& a,const FourierSeries& b){
    validate(a);validate(b);if(a.angularFrequency!=b.angularFrequency)throw std::invalid_argument("Fourier products require the same declared fundamental frequency");
    if(a.coefficients.size() && b.coefficients.size()>10000000/a.coefficients.size())throw std::invalid_argument("Fourier convolution exceeds 10000000 term products");
    FourierSeries out{a.angularFrequency,{}};
    for(const auto& [n,c]:a.coefficients)for(const auto& [m,d]:b.coefficients)out.coefficients[n+m]+=c*d;
    validate(out);return out;
}
FourierSeries projectCyclicSector(const FourierSeries& s,unsigned degree,int residue){
    validate(s);if(!degree||degree>1000000)throw std::invalid_argument("Cyclic degree must be 1..1000000");
    const int d=int(degree),r=((residue%d)+d)%d;FourierSeries out{s.angularFrequency,{}};
    for(const auto& [n,c]:s.coefficients)if(((n%d)+d)%d==r)out.coefficients[n]=c;return out;
}
double meanSquare(const FourierSeries& s){validate(s);double out=0;for(const auto& item:s.coefficients)out+=std::norm(item.second);if(!std::isfinite(out))throw std::overflow_error("Parseval sum overflow");return out;}
double omittedUniformBound(const FourierSeries& s,int maximum){
    validate(s);if(maximum<0)throw std::invalid_argument("Nonnegative truncation required");double out=0;
    for(const auto& [n,c]:s.coefficients)if(std::abs(n)>maximum)out+=std::abs(c);
    if(!std::isfinite(out))throw std::overflow_error("Fourier tail sum overflow");return out;
}
FourierSeries periodicFourierSeries(const std::vector<std::complex<double>>& samples,double frequency,int maximum){
    if(samples.empty()||maximum<0||maximum>1000000||std::size_t(maximum)*2>=samples.size())throw std::invalid_argument("Need 2*maximumHarmonic < sample count");
    FourierSeries out{frequency,{}};validate(out);
    for(auto c:samples)if(!std::isfinite(c.real())||!std::isfinite(c.imag()))throw std::invalid_argument("Fourier samples must be finite");
    const auto transform=core::discreteFourierTransform(samples);
    for(int n=-maximum;n<=maximum;++n){const auto i=n<0?samples.size()-std::size_t(-n):std::size_t(n);out.coefficients[n]=transform[i]/double(samples.size());}
    validate(out);return out;
}
std::vector<std::complex<double>> uniformFourierSamples(const FourierSeries& s,std::size_t count){
    validate(s);if(count<3||count>1000000)throw std::invalid_argument("Fourier grid must contain 3..1000000 samples");
    std::vector<std::complex<double>> c(count);
    for(const auto& [n,value]:s.coefficients){
        if(2*std::size_t(std::abs(n))>=count)throw std::invalid_argument("Fourier grid cannot resolve supplied harmonics");
        const auto i=n<0?count-std::size_t(-n):std::size_t(n);c[i]=value*double(count);
    }return core::discreteFourierTransform(c,true);
}
std::vector<std::complex<double>> periodicGridValues(const FourierSeries&s,std::size_t count,double phase){
    validate(s);if(count<3||count>1000000||!std::isfinite(phase))throw std::invalid_argument("Invalid Fourier point grid");
    phase=std::remainder(phase,2*std::acos(-1.));std::vector<std::complex<double>> c(count);const int N=int(count);
    for(const auto&[n,value]:s.coefficients)c[std::size_t((n%N+N)%N)]+=value*std::exp(std::complex<double>(0,n*phase))*double(count);
    for(auto z:c)finite(z);auto out=core::discreteFourierTransform(c,true);for(auto z:out)finite(z);return out;
}
}
