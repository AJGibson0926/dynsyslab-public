#include <dynsyslab/math/KoopmanFrequencyResponse.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace dynsyslab::math {
namespace {
using C=std::complex<double>;
const double pi=std::acos(-1.);
bool finiteMatrix(const ComplexMatrix& m){
    return m.data.size()==m.rows*m.cols&&std::all_of(m.data.begin(),m.data.end(),[](C z){return std::isfinite(z.real())&&std::isfinite(z.imag());});
}
ComplexVector inputColumn(const BilinearLiftedModel& m,const ComplexVector& d){
    ComplexVector b(m.A.rows);
    for(std::size_t j=0;j<m.B.cols;++j)for(std::size_t i=0;i<m.A.rows;++i)b[i]+=m.B(i,j)*d[j];
    return b;
}
ComplexMatrix inputCoupling(const BilinearLiftedModel& m,const ComplexVector& d){
    ComplexMatrix out(m.A.rows,m.A.rows);
    for(std::size_t j=0;j<m.N.size();++j)for(std::size_t k=0;k<out.data.size();++k)out.data[k]+=d[j]*m.N[j].data[k];
    return out;
}
ComplexVector directionOrDefault(const BilinearLiftedModel& m,const ComplexVector& d){
    if(d.empty()){ComplexVector e(m.B.cols);e.at(0)=1;return e;}
    if(d.size()!=m.B.cols)throw std::invalid_argument("Input direction size mismatch");
    for(const auto& z:d)if(!std::isfinite(z.real())||!std::isfinite(z.imag()))throw std::invalid_argument("Input direction must be finite");
    return d;
}
void requireFrequency(double w){if(!(w>0)||!std::isfinite(w))throw std::invalid_argument("Frequency must be positive");}
double scaleOf(const ComplexMatrix& a){return std::max(1.,frobeniusNorm(a));}
}
void validateBilinearModel(const BilinearLiftedModel& m){
    const auto n=m.A.rows;
    if(n==0||m.A.cols!=n||!finiteMatrix(m.A))throw std::invalid_argument("A must be a finite nonempty square matrix");
    if(m.B.rows!=n||m.B.cols==0||!finiteMatrix(m.B))throw std::invalid_argument("B must be finite with one row per state");
    if(m.C.cols!=n||m.C.rows==0||!finiteMatrix(m.C))throw std::invalid_argument("C must be finite with one column per state");
    if(!m.N.empty()&&m.N.size()!=m.B.cols)throw std::invalid_argument("N needs one matrix per input");
    for(const auto& x:m.N)if(x.rows!=n||x.cols!=n||!finiteMatrix(x))throw std::invalid_argument("Each N must be finite and n-by-n");
}
HarmonicResponse harmonicResponse(const BilinearLiftedModel& m,double w,const HarmonicResponseOptions& o){
    validateBilinearModel(m);requireFrequency(w);
    if(o.maximumOrder<1||o.maximumOrder>64)throw std::invalid_argument("Harmonic order must be 1..64");
    if(!(o.singularityTolerance>0))throw std::invalid_argument("Singularity tolerance must be positive");
    const auto d=directionOrDefault(m,o.direction);const auto b=inputColumn(m,d);
    const auto coupling=m.N.empty()?ComplexMatrix(m.A.rows,m.A.rows):inputCoupling(m,d);
    HarmonicResponse out;out.frequency=w;
    const auto eigenvalues=complexEigenvalues(m.A);
    out.hurwitz=std::all_of(eigenvalues.begin(),eigenvalues.end(),[](C z){return z.real()<0;});
    if(!out.hurwitz)out.warnings.push_back("A is not Hurwitz: a steady periodic response need not exist or attract, so H_n is formal");
    bool resonant=false;const double scale=scaleOf(m.A);
    ComplexVector z=b;
    for(int n=1;n<=o.maximumOrder;++n){
        const C s(0,n*w);auto shifted=shiftedMatrix(m.A,s);
        const double sigma=complexSingularValues(shifted).back();
        out.sigmaMin.push_back(sigma);
        if(sigma<=o.singularityTolerance*scale*std::max(1.,n*w)){
            resonant=true;out.warnings.push_back("Order "+std::to_string(n)+": i n w is an eigenvalue of A (external resonance); H_n does not exist");
            break;
        }
        ComplexLu lu(std::move(shifted));
        ComplexVector rhs=n==1?b:product(coupling,z);
        z=lu.solve(std::move(rhs));
        out.states.push_back(z);out.outputs.push_back(product(m.C,z));
    }
    out.valid=out.hurwitz&&!resonant;
    if(out.states.size()>=2){
        const double last=euclideanNorm(out.states.back()),prior=euclideanNorm(out.states[out.states.size()-2]);
        if(last>0&&prior>0)out.seriesRadius=prior/last;
    }
    return out;
}
ComplexVector volterraKernel(const BilinearLiftedModel& m,const ComplexVector& s,const std::vector<std::size_t>& channels){
    validateBilinearModel(m);const auto order=s.size();
    if(order<1||order>8)throw std::invalid_argument("Volterra order must be 1..8");
    if(!channels.empty()&&channels.size()!=order)throw std::invalid_argument("Volterra channel list size mismatch");
    std::vector<std::size_t> channel(order,0);
    if(!channels.empty())channel=channels;
    for(auto c:channel)if(c>=m.B.cols)throw std::invalid_argument("Volterra channel out of range");
    if(order>1&&m.N.empty()){return ComplexVector(m.C.rows);}
    const std::size_t subsets=std::size_t(1)<<order;const auto n=m.A.rows;
    std::vector<ComplexVector> table(subsets);
    for(std::size_t mask=1;mask<subsets;++mask){
        C total{};for(std::size_t j=0;j<order;++j)if(mask>>j&1)total+=s[j];
        ComplexVector rhs(n);
        for(std::size_t j=0;j<order;++j){if(!(mask>>j&1))continue;
            const auto rest=mask^(std::size_t(1)<<j);
            if(rest==0){for(std::size_t i=0;i<n;++i)rhs[i]+=m.B(i,channel[j]);}
            else{const auto part=product(m.N[channel[j]],table[rest]);for(std::size_t i=0;i<n;++i)rhs[i]+=part[i];}
        }
        ComplexLu lu(shiftedMatrix(m.A,total));
        if(lu.singular())throw std::domain_error("Volterra kernel evaluated at a pole");
        table[mask]=lu.solve(std::move(rhs));
    }
    double factorial=1;for(std::size_t k=2;k<=order;++k)factorial*=double(k);
    auto out=product(m.C,table[subsets-1]);for(auto& x:out)x/=factorial;
    return out;
}
RealHarmonicBalance realHarmonicBalance(const BilinearLiftedModel& m,double w,double a,int K,const ComplexVector& direction){
    validateBilinearModel(m);requireFrequency(w);
    if(!std::isfinite(a))throw std::invalid_argument("Amplitude must be finite");
    if(K<1||K>256)throw std::invalid_argument("Harmonic bound must be 1..256");
    const auto d=directionOrDefault(m,direction);const auto b=inputColumn(m,d);
    const auto coupling=m.N.empty()?ComplexMatrix(m.A.rows,m.A.rows):inputCoupling(m,d);
    const auto n=m.A.rows,blocks=std::size_t(2*K+1),size=n*blocks;
    if(size>4096)throw std::invalid_argument("Harmonic balance system is too large");
    ComplexMatrix system(size,size);ComplexVector rhs(size);
    const auto block=[&](int k){return std::size_t(k+K)*n;};
    for(int k=-K;k<=K;++k){const auto at=block(k);
        for(std::size_t i=0;i<n;++i){
            for(std::size_t j=0;j<n;++j)system(at+i,at+j)=-m.A(i,j);
            system(at+i,at+i)+=C(0,k*w);
            if(k>-K)for(std::size_t j=0;j<n;++j)system(at+i,block(k-1)+j)=-.5*a*coupling(i,j);
            if(k<K)for(std::size_t j=0;j<n;++j)system(at+i,block(k+1)+j)=-.5*a*coupling(i,j);
            if(k==1||k==-1)rhs[at+i]=.5*a*b[i];
        }
    }
    ComplexLu lu(std::move(system));if(lu.singular())throw std::domain_error("Harmonic balance system is singular: resonance at this frequency");
    const auto x=lu.solve(std::move(rhs));
    RealHarmonicBalance out;out.maximumHarmonic=K;
    for(int k=-K;k<=K;++k){const auto at=block(k);
        ComplexVector z(x.begin()+std::ptrdiff_t(at),x.begin()+std::ptrdiff_t(at+n));
        out.outputs[k]=product(m.C,z);out.states[k]=std::move(z);}
    auto edge=product(coupling,out.states[K]);auto low=product(coupling,out.states[-K]);
    out.truncationResidual=.5*std::abs(a)*std::max(euclideanNorm(edge),euclideanNorm(low));
    return out;
}
BilinearFloquet bilinearFloquet(const BilinearLiftedModel& m,double w,double a,bool complexInput,const ComplexVector& direction,std::size_t steps){
    validateBilinearModel(m);requireFrequency(w);
    if(!std::isfinite(a))throw std::invalid_argument("Amplitude must be finite");
    if(steps<16||steps>1048576)throw std::invalid_argument("Floquet steps per period out of range");
    const auto d=directionOrDefault(m,direction);
    const auto coupling=m.N.empty()?ComplexMatrix(m.A.rows,m.A.rows):inputCoupling(m,d);
    const auto n=m.A.rows;const double T=2*pi/w,h=T/double(steps);
    const auto input=[&](double t){return complexInput?C(a)*std::exp(C(0,w*t)):C(a*std::cos(w*t));};
    const auto rate=[&](double t,const ComplexMatrix& x){
        ComplexMatrix g(n,n);const C u=input(t);
        for(std::size_t k=0;k<g.data.size();++k)g.data[k]=m.A.data[k]+u*coupling.data[k];
        return product(g,x);};
    const auto combine=[&](const ComplexMatrix& x,const ComplexMatrix& v,double f){
        auto out=x;for(std::size_t k=0;k<out.data.size();++k)out.data[k]+=f*v.data[k];return out;};
    auto phi=ComplexMatrix::identity(n);
    for(std::size_t k=0;k<steps;++k){const double t=double(k)*h;
        const auto k1=rate(t,phi),k2=rate(t+.5*h,combine(phi,k1,.5*h)),k3=rate(t+.5*h,combine(phi,k2,.5*h)),k4=rate(t+h,combine(phi,k3,h));
        for(std::size_t e=0;e<phi.data.size();++e)phi.data[e]+=h/6*(k1.data[e]+2.*k2.data[e]+2.*k3.data[e]+k4.data[e]);
    }
    BilinearFloquet out;out.multipliers=complexEigenvalues(phi);
    for(const auto& z:out.multipliers)out.spectralRadius=std::max(out.spectralRadius,std::abs(z));
    out.stable=out.spectralRadius<1;
    return out;
}
SimplePoleCheck simplePoleCheck(const LaplaceFunction& F,C pole,double radius,double tolerance){
    if(!(tolerance>0))throw std::invalid_argument("Tolerance must be positive");
    CauchyOptions o;o.maximumOrder=2;const auto e=cauchyLaurent(F,pole,radius,o);
    SimplePoleCheck out;out.residue=e.laurent[0];out.secondOrderCoefficient=e.laurent[1];out.converged=e.converged;
    out.simple=e.converged&&std::abs(e.laurent[1])<=tolerance*std::max(std::abs(e.laurent[0]),std::numeric_limits<double>::min())&&e.higherOrderRatio<=tolerance;
    return out;
}
}
