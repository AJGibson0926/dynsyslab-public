#include <dynsyslab/math/InverseScattering.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>

// Marchenko reflectionless determinant, same KdV convention as Parker & Page,
// https://arxiv.org/abs/1903.04247, equations (24)-(25); log(c^2) rather than c.
// Inverse scattering theory: Gardner/Greene/Kruskal/Miura,
// https://doi.org/10.1103/PhysRevLett.19.1095.
namespace dynsyslab::math {
namespace {
void check(const ReflectionlessData&d){const auto n=d.kappa.size();
    if(n>64||d.logNormingSquared.size()!=n)throw std::invalid_argument("Reflectionless data need 0..64 paired poles and norming constants");
    for(std::size_t i=0;i<n;++i){if(!std::isfinite(d.kappa[i])||d.kappa[i]<=0||!std::isfinite(d.logNormingSquared[i]))throw std::invalid_argument("Positive finite kappa and finite log norming required");
        for(std::size_t j=0;j<i;++j)if(d.kappa[i]==d.kappa[j])throw std::invalid_argument("Distinct simple bound-state poles required");}}
}
ReflectionlessData evolveKdV(const ReflectionlessData&d,double t){check(d);if(!std::isfinite(t))throw std::invalid_argument("Finite time required");auto r=d;
    for(std::size_t i=0;i<r.kappa.size();++i){r.logNormingSquared[i]+=8*std::pow(r.kappa[i],3)*t;if(!std::isfinite(r.logNormingSquared[i]))throw std::overflow_error("Norming evolution overflow");}return r;
}
std::complex<double> evolveKdVReflection(std::complex<double>r,double k,double t){const double phase=8*k*k*k*t;
    if(!std::isfinite(phase)||!std::isfinite(r.real())||!std::isfinite(r.imag()))throw std::invalid_argument("Finite scattering coefficient, wavenumber and time required");return r*std::exp(std::complex<double>(0,phase));
}
double reflectionlessPotential(const ReflectionlessData&d,double x){check(d);if(!std::isfinite(x))throw std::invalid_argument("Finite position required");const auto n=d.kappa.size();if(!n)return 0;
    if(n==1){const double z=std::abs(d.kappa[0]*x-(d.logNormingSquared[0]-std::log(2*d.kappa[0]))/2);
        const double a=std::exp(-2*z);return -8*d.kappa[0]*d.kappa[0]*a/std::pow(1+a,2);}
    std::vector<double>M(n*n),D(n*n),DD(n*n),L(n*n);
    for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j){const double k=d.kappa[i]+d.kappa[j],logA=(d.logNormingSquared[i]+d.logNormingSquared[j])/2-k*x-std::log(k);
        if(logA>600)throw std::overflow_error("Reflectionless determinant requires a rescaled chart at this position/time");
        const double a=std::exp(logA);M[i*n+j]=a+(i==j?1.:0.);D[i*n+j]=-k*a;DD[i*n+j]=k*k*a;}
    for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<=i;++j){double v=M[i*n+j];for(std::size_t k=0;k<j;++k)v-=L[i*n+k]*L[j*n+k];
        if(i==j){if(v<=0||!std::isfinite(v))throw std::runtime_error("Reflectionless determinant lost positive definiteness; restrict/rescale the chart");L[i*n+j]=std::sqrt(v);}else L[i*n+j]=v/L[j*n+j];}
    auto solve=[&](const std::vector<double>&B){std::vector<double>Y(n*n),X(n*n);
        for(std::size_t c=0;c<n;++c){for(std::size_t i=0;i<n;++i){double v=B[i*n+c];for(std::size_t j=0;j<i;++j)v-=L[i*n+j]*Y[j*n+c];Y[i*n+c]=v/L[i*n+i];}
            for(std::size_t i=n;i-->0;){double v=Y[i*n+c];for(std::size_t j=i+1;j<n;++j)v-=L[j*n+i]*X[j*n+c];X[i*n+c]=v/L[i*n+i];}}return X;};
    const auto A=solve(D),B=solve(DD);double derivative=0;
    for(std::size_t i=0;i<n;++i){derivative+=B[i*n+i];for(std::size_t j=0;j<n;++j)derivative-=A[i*n+j]*A[j*n+i];}
    if(!std::isfinite(derivative))throw std::overflow_error("Reflectionless reconstruction overflow");return -2*derivative;
}
}
