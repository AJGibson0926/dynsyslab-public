#include <dynsyslab/math/PronyIdentification.hpp>
#include <dynsyslab/math/ResolventLeastSquares.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace dynsyslab::math {
namespace {
using C=std::complex<double>;
const double pi=std::acos(-1.);
// Eigenvalues of pinv(W1) W2 from the leading M right singular vectors. The row space of the
// Hankel matrix is spanned by conj(V), and the Vandermonde shift gives W2=W1 T^{-1} diag(mu) T.
ComplexVector pencilPoles(const ComplexMatrix& V,std::size_t M){
    const auto rows=V.rows;ComplexMatrix W1(rows-1,M),W2(rows-1,M);
    for(std::size_t i=0;i+1<rows;++i)for(std::size_t j=0;j<M;++j){W1(i,j)=std::conj(V(i,j));W2(i,j)=std::conj(V(i+1,j));}
    return complexEigenvalues(leastSquares(W1,W2,1e-14).x);
}
struct Fit {ComplexVector residues;double relative{},condition{};};
Fit fitResidues(const ComplexVector& y,const ComplexVector& mu){
    ComplexMatrix A(y.size(),mu.size()),b(y.size(),1);
    for(std::size_t j=0;j<mu.size();++j){C p=1;for(std::size_t k=0;k<y.size();++k){A(k,j)=p;p*=mu[j];}}
    for(std::size_t k=0;k<y.size();++k)b(k,0)=y[k];
    const auto s=leastSquares(A,b,1e-14);Fit f;f.residues.resize(mu.size());f.condition=s.condition;
    for(std::size_t j=0;j<mu.size();++j)f.residues[j]=s.x(j,0);
    const auto r=product(A,f.residues);double e=0,t=0;for(std::size_t k=0;k<y.size();++k){e+=std::norm(r[k]-y[k]);t+=std::norm(y[k]);}
    f.relative=t>0?std::sqrt(e/t):std::sqrt(e);
    return f;
}
}
PronyResult pronyIdentify(const ComplexVector& y,const PronyOptions& o){
    const auto N=y.size();if(N<4||N>200000)throw std::invalid_argument("Prony needs 4..200000 samples");
    for(const auto& z:y)if(!std::isfinite(z.real())||!std::isfinite(z.imag()))throw std::invalid_argument("Samples must be finite");
    if(!(o.sampleInterval>0)||!std::isfinite(o.sampleInterval))throw std::invalid_argument("Positive sample interval required");
    if(!(o.rankTolerance>0)||!(o.clusterFactor>=0)||!(o.gapThreshold>1))throw std::invalid_argument("Invalid Prony tolerances");
    if(std::all_of(y.begin(),y.end(),[](C z){return z==C{};}))throw std::invalid_argument("Prony needs a nonzero signal");
    const std::size_t L=o.pencilParameter?o.pencilParameter:std::min<std::size_t>(N/2,128);
    if(L<1||L>512||L+1>=N)throw std::invalid_argument("Pencil parameter must satisfy 1<=L<=512 and L+1<N");
    ComplexMatrix H(N-L,L+1);
    for(std::size_t i=0;i<N-L;++i)for(std::size_t j=0;j<=L;++j)H(i,j)=y[i+j];
    const auto svd=thinSvd(H);const auto& s=svd.singularValues;const auto limit=std::min(L,N-L);
    PronyResult out;out.singularValues=s;out.pencilParameter=L;
    std::size_t M=0;
    if(o.order){M=*o.order;out.orderSource="declared";if(M<1||M>limit)throw std::invalid_argument("Declared order must lie in 1..min(L,N-L)");}
    else{out.orderSource="rank";for(double v:s)if(v>o.rankTolerance*s[0])++M;M=std::min(std::max<std::size_t>(M,1),limit);}
    out.order=M;
    const auto mu=pencilPoles(svd.V,M);auto fit=fitResidues(y,mu);
    const std::size_t scanTop=std::min<std::size_t>(limit,std::min<std::size_t>(M+3,16));
    for(std::size_t m=1;m<=scanTop;++m){
        if(m==M){out.orderScan.push_back({m,fit.relative});continue;}
        try{out.orderScan.push_back({m,fitResidues(y,pencilPoles(svd.V,m)).relative});}catch(const std::exception&){}
    }
    out.relativeResidual=fit.relative;out.vandermondeCondition=fit.condition;
    out.singularGap=M<s.size()?(s[M]>0?s[M-1]/s[M]:std::numeric_limits<double>::infinity()):0;
    std::vector<std::size_t> order(M);for(std::size_t j=0;j<M;++j)order[j]=j;
    std::stable_sort(order.begin(),order.end(),[&](std::size_t a,std::size_t b){return std::abs(fit.residues[a])>std::abs(fit.residues[b]);});
    out.model.sampleInterval=o.sampleInterval;std::size_t dropped=0;
    for(auto j:order){
        if(std::abs(mu[j])<1e-12){++dropped;continue;}
        out.poles.push_back({mu[j],fit.residues[j],C(std::log(std::abs(mu[j])),std::arg(mu[j]))/o.sampleInterval});
        out.model.poles.push_back({mu[j],{fit.residues[j]}});
    }
    if(dropped)out.warnings.push_back(std::to_string(dropped)+" pole(s) at z=0 dropped: they have no continuous-time eigenvalue");
    // Separation in log(mu), modulo 2 pi on the angle, against the Rayleigh limit 2 pi/N.
    out.rayleighLimit=2*pi/double(N);out.minimumSeparation=std::numeric_limits<double>::infinity();
    for(std::size_t a=0;a<out.poles.size();++a)for(std::size_t b=a+1;b<out.poles.size();++b){
        const C d=std::log(out.poles[a].multiplier)-std::log(out.poles[b].multiplier);
        out.minimumSeparation=std::min(out.minimumSeparation,std::hypot(d.real(),std::remainder(d.imag(),2*pi)));
    }
    out.clustered=out.poles.size()>1&&out.minimumSeparation<o.clusterFactor*out.rayleighLimit;
    if(out.clustered)out.warnings.push_back("Poles closer than the record's Rayleigh limit 2*pi/N: individual poles and residues are unreliable; this may be a repeated pole or a band");
    const bool noGap=M>=limit||(M<s.size()&&s[M]>0&&out.singularGap<o.gapThreshold);
    out.nearContinuous=noGap&&fit.relative>1e-6;
    if(M>=limit)out.warnings.push_back("Order reached the pencil limit: the signal is not low-rank at this tolerance");
    else if(noGap)out.warnings.push_back("No clear singular-value gap at the selected order");
    if(out.nearContinuous)out.warnings.push_back("Signal looks near-continuous in spectrum: finite poles only approximate it, and do not identify Koopman eigenvalues");
    if(fit.condition>1e10)out.warnings.push_back("Ill-conditioned Vandermonde residue fit");
    for(const auto& p:out.poles)if(std::abs(p.multiplier)>1+1e-6){out.warnings.push_back("Pole outside the unit circle: growing mode or a noise artifact");break;}
    out.warnings.push_back("Eigenvalues are on the principal branch: frequencies above pi/h are aliased and cannot be recovered from these samples");
    return out;
}
}
