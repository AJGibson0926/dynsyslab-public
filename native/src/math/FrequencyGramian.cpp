#include <dynsyslab/math/FrequencyGramian.hpp>
#ifdef DYNSYSLAB_MATH_NUMERICS
#include <dynsyslab/core/StateSpace.hpp>
#include <dynsyslab/core/SystemNorms.hpp>
#endif
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
namespace dynsyslab::math {
namespace {
using Cx=std::complex<double>;
const double pi=std::acos(-1.);
// Eight-point Gauss-Legendre rule on [-1,1] (Abramowitz & Stegun, Table 25.4).
const std::array<double,4> gaussX{.1834346424956498,.5255324099163290,.7966664774136267,.9602898564975363};
const std::array<double,4> gaussW{.3626837833783620,.3137066458778873,.2223810344533745,.1012285362903763};
void validate(const ComplexMatrix& A,const ComplexMatrix& B,const ComplexMatrix& C){
    if(A.rows==0||A.rows!=A.cols||A.rows>256)throw std::invalid_argument("A must be square with dimension 1..256");
    if(B.rows!=A.rows||B.cols==0)throw std::invalid_argument("B needs one row per state");
    if(C.cols!=A.rows||C.rows==0)throw std::invalid_argument("C needs one column per state");
    for(const auto* m:{&A,&B,&C})for(const auto& z:m->data)if(!std::isfinite(z.real())||!std::isfinite(z.imag()))throw std::invalid_argument("Matrices must be finite");
}
void addScaled(ComplexMatrix& target,const ComplexMatrix& x,double weight){for(std::size_t k=0;k<x.data.size();++k)target.data[k]+=weight*x.data[k];}
}
FrequencyGramianResult frequencyGramians(const ComplexMatrix& A,const ComplexMatrix& B,const ComplexMatrix& C,const FrequencyGramianOptions& o){
    validate(A,B,C);
    if(!(o.tolerance>0)||o.maximumLevels<1||o.maximumLevels>12)throw std::invalid_argument("Invalid Gramian quadrature options");
    const auto eigenvalues=complexEigenvalues(A);
    std::vector<double> mags;
    for(const auto& z:eigenvalues){
        if(!(z.real()<0))throw std::domain_error("Frequency-domain Gramians need a Hurwitz A");
        mags.push_back(std::abs(z));
    }
    std::nth_element(mags.begin(),mags.begin()+std::ptrdiff_t(mags.size()/2),mags.end());
    const double c=std::max(mags[mags.size()/2],1e-12);
    // Breakpoints in phi=atan(w/c), graded geometrically at each pole frequency.
    std::vector<double> breaks{-pi/2,pi/2};
    for(const auto& z:eigenvalues){const double w=z.imag(),d=std::max(-z.real(),1e-14*c);
        breaks.push_back(std::atan(w/c));
        for(double f=1;f<1e6*c/d&&f<1e9;f*=2){breaks.push_back(std::atan((w+f*d)/c));breaks.push_back(std::atan((w-f*d)/c));}
    }
    std::sort(breaks.begin(),breaks.end());
    breaks.erase(std::unique(breaks.begin(),breaks.end(),[](double a,double b){return std::abs(a-b)<1e-13;}),breaks.end());
    const auto BB=product(B,adjoint(B)),CC=product(adjoint(C),C);
    FrequencyGramianResult out;
    const auto integrate=[&](std::size_t subdivisions,ComplexMatrix& wc,ComplexMatrix& wo){
        wc=ComplexMatrix(A.rows,A.rows);wo=ComplexMatrix(A.rows,A.rows);
        for(std::size_t b=0;b+1<breaks.size();++b){const double width=(breaks[b+1]-breaks[b])/double(subdivisions);
            for(std::size_t p=0;p<subdivisions;++p){const double mid=breaks[b]+(double(p)+.5)*width,r=.5*width;
                for(int side=-1;side<=1;side+=2)for(std::size_t a=0;a<4;++a){
                    const double phi=mid+side*r*gaussX[a],w=c*std::tan(phi),cosine=std::cos(phi);
                    const double jacobian=c/(cosine*cosine),weight=r*gaussW[a]*jacobian/(2*pi);
                    const ComplexLu lu(shiftedMatrix(A,Cx(0,w)));
                    // R BB^H R^H = R (R BB^H)^H, and R^H CC R = R^H (R^H CC)^H.
                    const auto X=lu.solve(BB);addScaled(wc,lu.solve(adjoint(X)),weight);
                    const auto Y=lu.solveAdjoint(CC);addScaled(wo,lu.solveAdjoint(adjoint(Y)),weight);
                    ++out.evaluations;
                }
            }
        }
    };
    ComplexMatrix wc,wo,nextC,nextO;integrate(1,wc,wo);std::size_t subdivisions=1;double change=1;
    for(std::size_t level=1;level<o.maximumLevels;++level){
        subdivisions*=2;integrate(subdivisions,nextC,nextO);
        ComplexMatrix dc=nextC,dO=nextO;addScaled(dc,wc,-1);addScaled(dO,wo,-1);
        change=std::max(frobeniusNorm(dc)/std::max(frobeniusNorm(nextC),1e-300),frobeniusNorm(dO)/std::max(frobeniusNorm(nextO),1e-300));
        wc=nextC;wo=nextO;
        if(change<=o.tolerance){out.converged=true;break;}
    }
    // The quadrature is Hermitian only to rounding: symmetrize.
    for(auto* m:{&wc,&wo})for(std::size_t i=0;i<m->rows;++i)for(std::size_t j=i;j<m->cols;++j){
        const Cx v=.5*((*m)(i,j)+std::conj((*m)(j,i)));(*m)(i,j)=v;(*m)(j,i)=std::conj(v);}
    const auto outputs=product(product(C,wc),adjoint(C));double trace=0;
    for(std::size_t i=0;i<outputs.rows;++i)trace+=outputs(i,i).real();
    out.h2Norm=std::sqrt(std::max(trace,0.));out.relativeChange=change;
    out.controllability=std::move(wc);out.observability=std::move(wo);
    return out;
}
ResolventPeak resolventHInfinity(const ComplexMatrix& A,const ComplexMatrix& B,const ComplexMatrix& C){
    validate(A,B,C);
#ifdef DYNSYSLAB_MATH_NUMERICS
    const auto realify=[](const ComplexMatrix& m){
        std::vector<double> out(4*m.rows*m.cols);const auto R=2*m.rows,K=2*m.cols;
        for(std::size_t i=0;i<m.rows;++i)for(std::size_t j=0;j<m.cols;++j){const auto z=m(i,j);
            out[i*K+j]=z.real();out[i*K+j+m.cols]=-z.imag();out[(i+m.rows)*K+j]=z.imag();out[(i+m.rows)*K+j+m.cols]=z.real();}
        (void)R;return out;};
    const auto model=core::makeStateSpace(realify(A),realify(B),realify(C),{},2*A.rows,2*B.cols,2*C.rows);
    const auto norm=core::hInfinityNorm(model);
    const auto gain=[&](double w){const ComplexLu lu(shiftedMatrix(A,Cx(0,w)));
        return complexSingularValues(product(C,lu.solve(B))).front();};
    ResolventPeak out;out.value=norm.value;out.iterations=norm.iterations;
    out.frequency=std::isfinite(norm.value)?(gain(norm.peakFrequency)>=gain(-norm.peakFrequency)?norm.peakFrequency:-norm.peakFrequency):0;
    return out;
#else
    throw std::runtime_error("H-infinity peak needs the optional numerics library");
#endif
}
}
