#include <dynsyslab/math/TorusFourier.hpp>
#include <dynsyslab/core/DiscreteFourier.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>

// Fourier analysis on T^d; Arnold, Mathematical Methods of Classical Mechanics,
// chapters 10 and appendices 7-8, https://doi.org/10.1007/978-1-4757-2063-1.
// Cohomological equation used by parameterization methods:
// Figueras/Haro/Luque https://arxiv.org/abs/1601.00084.
namespace dynsyslab::math {
namespace {
bool finite(std::complex<double> z){return std::isfinite(z.real())&&std::isfinite(z.imag());}
void check(const TorusFourierSeries&s){
    if(s.angularFrequencies.empty()||s.angularFrequencies.size()>8||s.coefficients.size()>1000000)throw std::invalid_argument("Invalid torus size");
    for(auto w:s.angularFrequencies)if(!std::isfinite(w))throw std::invalid_argument("Finite frequencies required");
    for(const auto&[k,c]:s.coefficients){if(k.size()!=s.angularFrequencies.size()||!finite(c))throw std::invalid_argument("Invalid torus coefficient");
        for(auto n:k)if(n<-1000000||n>1000000)throw std::invalid_argument("Harmonic index too large");}
}
}
double harmonicFrequency(const FourierIndex&k,const std::vector<double>&w){
    if(k.size()!=w.size())throw std::invalid_argument("Frequency dimension mismatch");double v=0;
    for(std::size_t j=0;j<k.size();++j)v+=k[j]*w[j];if(!std::isfinite(v))throw std::invalid_argument("Nonfinite harmonic frequency");return v;
}
TorusFourierSeries torusFourierSeries(const std::vector<std::complex<double>>& samples,
    const std::vector<std::size_t>&shape,const std::vector<double>&w,const std::vector<int>&maximum){
    TorusFourierSeries out{w,{}};check(out);const auto d=w.size();
    if(shape.size()!=d||maximum.size()!=d)throw std::invalid_argument("Torus mesh dimension mismatch");
    std::size_t total=1;for(std::size_t j=0;j<d;++j){
        if(shape[j]<3||shape[j]>1000000||total>1000000/shape[j]||maximum[j]<0||2*std::size_t(maximum[j])>=shape[j])throw std::invalid_argument("Unresolved/oversized torus mesh");
        total*=shape[j];}
    if(samples.size()!=total)throw std::invalid_argument("Torus sample count mismatch");
    for(auto z:samples)if(!finite(z))throw std::invalid_argument("Nonfinite torus sample");auto c=samples;
    std::size_t stride=1;
    for(std::size_t axis=d;axis-->0;){const auto n=shape[axis],block=n*stride;
        for(std::size_t offset=0;offset<total;offset+=block)for(std::size_t tail=0;tail<stride;++tail){
            std::vector<std::complex<double>> line(n);for(std::size_t j=0;j<n;++j)line[j]=c[offset+tail+j*stride];
            auto transformed=core::discreteFourierTransform(line);for(std::size_t j=0;j<n;++j)c[offset+tail+j*stride]=transformed[j];}
        stride=block;}
    for(std::size_t flat=0;flat<total;++flat){std::size_t i=flat;FourierIndex k(d);bool keep=true;
        for(std::size_t axis=d;axis-->0;){const auto index=i%shape[axis];i/=shape[axis];
            k[axis]=index<=shape[axis]/2?int(index):int(index)-int(shape[axis]);keep=keep&&std::abs(k[axis])<=maximum[axis];}
        if(keep)out.coefficients.emplace(std::move(k),c[flat]/double(total));}
    check(out);return out;
}
std::complex<double> evaluate(const TorusFourierSeries&s,const std::vector<double>&phi){
    check(s);if(phi.size()!=s.angularFrequencies.size())throw std::invalid_argument("Phase dimension mismatch");
    std::vector<double> reduced;for(auto p:phi){if(!std::isfinite(p))throw std::invalid_argument("Finite phases required");reduced.push_back(std::remainder(p,2*std::acos(-1.)));}
    std::complex<double> out{};for(const auto&[k,c]:s.coefficients)out+=c*std::exp(std::complex<double>(0,harmonicFrequency(k,reduced)));
    if(!finite(out))throw std::overflow_error("Torus evaluation overflow");return out;
}
TorusFourierSeries torusTimeDerivative(const TorusFourierSeries&s){check(s);auto r=s;
    for(auto&[k,c]:r.coefficients)c*=std::complex<double>(0,harmonicFrequency(k,s.angularFrequencies));check(r);return r;}
TorusFourierSeries torusPartialDerivative(const TorusFourierSeries&s,std::size_t axis){check(s);
    if(axis>=s.angularFrequencies.size())throw std::invalid_argument("Torus derivative axis out of range");auto r=s;
    for(auto&[k,c]:r.coefficients)c*=std::complex<double>(0,k[axis]);check(r);return r;}
double meanSquare(const TorusFourierSeries&s){check(s);double v=0;for(const auto&[k,c]:s.coefficients)v+=std::norm(c);
    if(!std::isfinite(v))throw std::overflow_error("Torus norm overflow");return v;}
CohomologicalSolution solveCohomological(const TorusFourierSeries&g,double tol){check(g);
    if(!std::isfinite(tol)||tol<0)throw std::invalid_argument("Nonnegative finite divisor tolerance required");
    CohomologicalSolution r{{g.angularFrequencies,{}},{g.angularFrequencies,{}},0};
    for(const auto&[k,c]:g.coefficients){const double f=harmonicFrequency(k,g.angularFrequencies);
        if(std::abs(f)<=tol)r.retained.coefficients[k]=c;
        else{r.solution.coefficients[k]=c/std::complex<double>(0,f);r.maximumInverseDivisor=std::max(r.maximumInverseDivisor,1/std::abs(f));}}
    check(r.solution);return r;
}
}
