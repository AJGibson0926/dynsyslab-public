#include <dynsyslab/math/ForcedResponseEstimation.hpp>
#include <dynsyslab/math/PhaseCoverFourier.hpp>
#include <dynsyslab/math/ResolventLeastSquares.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace dynsyslab::math {
namespace {
using C=std::complex<double>;
const double pi=std::acos(-1.);
C monomial(const ComplexVector& x,C u,const std::vector<int>& e){
    C v=1;for(std::size_t i=0;i<x.size();++i)for(int p=0;p<e[i];++p)v*=x[i];
    for(int p=0;p<e[x.size()];++p)v*=u;
    return v;
}
void enumerate(std::size_t vars,int degree,std::vector<std::vector<int>>& out){
    std::vector<int> e(vars,0);
    const auto walk=[&](const auto& self,std::size_t i,int left)->void{
        if(i==vars){if(std::any_of(e.begin(),e.end(),[](int v){return v>0;}))out.push_back(e);return;}
        for(int v=0;v<=left;++v){e[i]=v;self(self,i+1,left-v);}
        e[i]=0;};
    walk(walk,0,degree);
}
std::size_t findMonomial(const std::vector<std::vector<int>>& d,const std::vector<int>& e){
    const auto it=std::find(d.begin(),d.end(),e);return it==d.end()?d.size():std::size_t(it-d.begin());
}
// (p,D) with n=p^D, or {0,0}.
std::pair<std::uint32_t,std::size_t> primePower(int n){
    for(int p=2;p*p<=n;++p)if(n%p==0){int m=n;std::size_t D=0;while(m%p==0){m/=p;++D;}
        return m==1?std::pair<std::uint32_t,std::size_t>{std::uint32_t(p),D}:std::pair<std::uint32_t,std::size_t>{0,0};}
    return {std::uint32_t(n),1};
}
}
ForcedHarmonicEstimate estimateForcedHarmonics(const SampledSignal& g,double w,double u0,int K,bool complexInput,double fraction){
    if(!(u0>0)||!std::isfinite(u0))throw std::invalid_argument("Forcing amplitude must be positive");
    ForcedHarmonicEstimate out;out.coefficients=estimateHarmonics(g,w,K,fraction);
    for(int n=1;n<=K;++n)if(complexInput)out.response[n]=out.coefficients.at(n)/std::pow(u0,n);
    if(!complexInput)out.warnings.push_back("Real input: c_n mixes Volterra orders, so H_n is not reported; use realHarmonicBalance or volterraKernel");
    const double start=g.times.front()+fraction*(g.times.back()-g.times.front());double e=0,t=0;
    for(std::size_t k=0;k<g.times.size();++k){if(g.times[k]<start)continue;C fit{};
        for(const auto& [n,c]:out.coefficients)fit+=c*std::exp(C(0,n*w*g.times[k]));
        e+=std::norm(g.values[k]-fit);t+=std::norm(g.values[k]);}
    out.tailRelativeResidual=t>0?std::sqrt(e/t):std::sqrt(e);
    if(out.tailRelativeResidual>1e-3)out.warnings.push_back("Tail is not explained by harmonics of w: not settled, subharmonic, quasi-periodic, or under-resolved");
    return out;
}
SkewDmdResult skewProductDmd(const SkewDmdInput& q){
    const auto T=q.states.size();if(T<4||q.inputs.size()!=T)throw std::invalid_argument("Skew DMD needs at least four matching snapshots");
    const auto d=q.states[0].size();if(d==0||q.outputState>=d)throw std::invalid_argument("Invalid state dimension or output state");
    for(const auto& x:q.states){if(x.size()!=d)throw std::invalid_argument("Snapshot size mismatch");
        for(const auto& z:x)if(!std::isfinite(z.real())||!std::isfinite(z.imag()))throw std::invalid_argument("Snapshots must be finite");}
    if(!(q.sampleInterval>0)||!(q.frequency>0)||!std::isfinite(q.sampleInterval)||!std::isfinite(q.frequency))throw std::invalid_argument("Positive sample interval and frequency required");
    if(q.maximumHarmonic<1||q.maximumHarmonic>16||q.maximumDegree<1||q.maximumDegree>8)throw std::invalid_argument("Harmonic must be 1..16 and degree 1..8");
    auto dict=q.dictionary;if(dict.empty())enumerate(d+1,q.maximumDegree,dict);
    for(const auto& e:dict){if(e.size()!=d+1)throw std::invalid_argument("Dictionary exponents need one entry per state plus u");
        if(std::any_of(e.begin(),e.end(),[](int v){return v<0||v>16;}))throw std::invalid_argument("Dictionary exponents must be 0..16");}
    const auto M=dict.size();if(M>256)throw std::invalid_argument("Dictionary exceeds 256 functions");
    SkewDmdResult out;out.dictionarySize=M;out.snapshots=T-1;
    if(T-1<M){out.warnings.push_back("Fewer snapshot pairs than dictionary functions");return out;}
    ComplexMatrix psi(T-1,M),next(T-1,M);std::vector<double> scale(M,0);
    for(std::size_t k=0;k<T;++k)for(std::size_t j=0;j<M;++j){const C v=monomial(q.states[k],q.inputs[k],dict[j]);
        if(k+1<T)psi(k,j)=v;if(k>0)next(k-1,j)=v;scale[j]=std::max(scale[j],std::abs(v));}
    for(auto& s:scale)if(!(s>0))s=1;
    // Equilibrated columns: Psi D^{-1} (D K^T D^{-1}) = Psi+ D^{-1}, rows are snapshots.
    ComplexMatrix P(T-1,M),Q(T-1,M);
    for(std::size_t k=0;k<T-1;++k)for(std::size_t j=0;j<M;++j){P(k,j)=psi(k,j)/scale[j];Q(k,j)=next(k,j)/scale[j];}
    const auto fit=leastSquares(P,Q,q.rankTolerance);out.rank=fit.rank;out.condition=fit.condition;
    ComplexMatrix K(M,M);for(std::size_t i=0;i<M;++i)for(std::size_t j=0;j<M;++j)K(i,j)=fit.x(j,i)*scale[i]/scale[j];
    double e=0,t=0;const auto predicted=product(P,fit.x);
    for(std::size_t k=0;k<Q.data.size();++k){e+=std::norm(predicted.data[k]-Q.data[k]);t+=std::norm(Q.data[k]);}
    out.closureResidual=t>0?std::sqrt(e/t):0;
    if(fit.rank<M){out.warnings.push_back("Snapshot matrix is rank deficient: steady-state-only data cannot identify the lift; include the transient or a smaller dictionary");return out;}
    std::vector<int> outputExponent(d+1,0);outputExponent[q.outputState]=1;const auto outputIndex=findMonomial(dict,outputExponent);
    if(outputIndex==M){out.warnings.push_back("Dictionary lacks the output coordinate");return out;}
    out.feasible=true;const auto spectrum=complexEigenvalues(K);
    for(int n=1;n<=q.maximumHarmonic;++n){
        SkewDmdHarmonic h;h.harmonic=n;const C mu=std::exp(C(0,n*q.frequency*q.sampleInterval));
        h.spectrumDistance=std::numeric_limits<double>::infinity();for(const auto& z:spectrum)h.spectrumDistance=std::min(h.spectrumDistance,std::abs(z-mu));
        std::vector<int> ue(d+1,0);ue[d]=n;const auto index=findMonomial(dict,ue);
        if(index==M){out.warnings.push_back("Dictionary lacks u^"+std::to_string(n));out.harmonics.push_back(h);continue;}
        ComplexMatrix A(M-1,M-1);ComplexVector b(M-1);
        const auto reduce=[&](std::size_t r){return r<index?r:r-1;};
        for(std::size_t r=0;r<M;++r){if(r==index)continue;b[reduce(r)]=K(r,index);
            for(std::size_t c=0;c<M;++c){if(c==index)continue;A(reduce(r),reduce(c))=(r==c?mu:C{})-K(r,c);}}
        ComplexVector v(M);v[index]=1;
        if(M>1){ComplexLu lu(std::move(A));if(lu.singular()){out.warnings.push_back("Mode solve singular at harmonic "+std::to_string(n));out.harmonics.push_back(h);continue;}
            const auto x=lu.solve(std::move(b));for(std::size_t r=0;r<M;++r)if(r!=index)v[r]=x[reduce(r)];}
        C top{};for(std::size_t c=0;c<M;++c)top+=K(index,c)*v[c];
        h.modeResidual=std::abs(mu-top);h.response=v[outputIndex];
        if(h.modeResidual>1e-6)out.warnings.push_back("Harmonic "+std::to_string(n)+": u^n is not an accurate eigenfunction of the fitted K");
        out.harmonics.push_back(h);
    }
    if(out.closureResidual>1e-6)out.warnings.push_back("Dictionary is not closed under the flow: H_n carries the closure error");
    return out;
}
SubharmonicEstimate subharmonicResponse(const SampledSignal& g,double w,double u0,int n,int maximumIndex,double fraction,double tolerance){
    if(n<2||n>4096)throw std::invalid_argument("Cover degree must be 2..4096");
    if(!(w>0)||!std::isfinite(w)||!(u0>0)||!std::isfinite(u0))throw std::invalid_argument("Positive frequency and amplitude required");
    if(maximumIndex<1||maximumIndex>4096||!(tolerance>0))throw std::invalid_argument("Invalid subharmonic bounds");
    SubharmonicEstimate out;out.cover=n;out.coefficients=estimateHarmonics(g,w/n,maximumIndex,fraction);
    double sub=0,total=0;
    for(const auto& [k,c]:out.coefficients){if(k==0)continue;const double e=std::norm(c);total+=e;if(k%n!=0)sub+=e;}
    out.subharmonicFraction=total>0?sub/total:0;out.detected=out.subharmonicFraction>tolerance;
    if(!out.detected)out.warnings.push_back("No subharmonic energy above tolerance: H_{1/n} is not reported");
    else{
        const C H=out.coefficients.at(1)/std::pow(u0,1./n);
        if(std::abs(H)>0){
            out.magnitude=std::abs(H);const double width=2*pi/n;double a=std::fmod(std::arg(H),width);if(a<0)a+=width;
            out.phaseModulo=a;out.invariantPower=std::pow(H,n);
            for(int j=0;j<n;++j)out.branches.push_back(H*std::exp(C(0,2*pi*j/n)));
        }else out.warnings.push_back("Subharmonic energy is not at k=1: H_{1/n} is zero; inspect the coefficients");
        out.warnings.push_back("H_{1/n} is defined only up to an n-th root of unity; the branch depends on the initial phase");
    }
    const auto [p,D]=primePower(n);
    if(p){
        std::map<FourierIndex,C> modes;for(const auto& [k,c]:out.coefficients)if(c!=C{})modes[FourierIndex{k}]=c;
        if(!modes.empty()){const auto r=phaseCoverForcing(p,D,w,{},modes,1e-12*w);out.coverLabelled=true;
            for(const auto& m:r.modes)out.labels.push_back({m.index[0],m.frequencyRatio,m.minimumDepth});}
    }else out.warnings.push_back("Cover degree is not a prime power: PhaseCoverFourier labels are unavailable");
    return out;
}
}
