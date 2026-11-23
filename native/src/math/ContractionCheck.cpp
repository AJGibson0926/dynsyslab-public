#include <dynsyslab/math/ContractionCheck.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace dynsyslab::math {
namespace {
using C=std::complex<double>;
std::vector<double> weightsOrUnit(const std::vector<double>& w,std::size_t n){
    if(w.empty())return std::vector<double>(n,1.);
    if(w.size()!=n)throw std::invalid_argument("One positive weight per state required");
    for(double x:w)if(!(x>0)||!std::isfinite(x))throw std::invalid_argument("Weights must be positive and finite");
    return w;
}
double largestHermitianEigenvalue(const ComplexMatrix& J){
    const auto n=J.rows;ComplexMatrix H(n,n);
    for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j)H(i,j)=.5*(J(i,j)+std::conj(J(j,i)));
    double top=-std::numeric_limits<double>::infinity();
    for(const auto& z:complexEigenvalues(H))top=std::max(top,z.real());
    return top;
}
std::pair<double,const char*> bestOf(const LogNorms& n){
    if(n.two<=n.one&&n.two<=n.infinity)return {n.two,"two"};
    return n.one<=n.infinity?std::pair<double,const char*>{n.one,"one"}:std::pair<double,const char*>{n.infinity,"infinity"};
}
}
LogNorms logarithmicNorms(const ComplexMatrix& J,const std::vector<double>& weights){
    if(J.rows==0||J.rows!=J.cols)throw std::invalid_argument("Square Jacobian required");
    for(const auto& z:J.data)if(!std::isfinite(z.real())||!std::isfinite(z.imag()))throw std::invalid_argument("Jacobian must be finite");
    const auto n=J.rows;const auto w=weightsOrUnit(weights,n);ComplexMatrix S(n,n);
    for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j)S(i,j)=J(i,j)*w[i]/w[j];
    LogNorms out;out.two=largestHermitianEigenvalue(S);out.one=out.infinity=-std::numeric_limits<double>::infinity();
    for(std::size_t i=0;i<n;++i){double row=S(i,i).real(),column=S(i,i).real();
        for(std::size_t j=0;j<n;++j)if(j!=i){row+=std::abs(S(i,j));column+=std::abs(S(j,i));}
        out.infinity=std::max(out.infinity,row);out.one=std::max(out.one,column);}
    return out;
}
SampledContractionReport sampledContraction(const SampledContractionInput& q){
    const auto D=q.lower.size();
    if(q.stateDimension<1||q.stateDimension>D||q.upper.size()!=D||!q.jacobian)throw std::invalid_argument("Contraction check needs a Jacobian and a box covering the state");
    for(std::size_t i=0;i<D;++i)if(!std::isfinite(q.lower[i])||!std::isfinite(q.upper[i])||q.lower[i]>q.upper[i])throw std::invalid_argument("Invalid box bounds");
    if(q.samplesPerAxis<2||q.samplesPerAxis>257)throw std::invalid_argument("Samples per axis must be 2..257");
    double total=1;for(std::size_t i=0;i<D;++i)total*=double(q.lower[i]==q.upper[i]?1:q.samplesPerAxis);
    if(total>2e6)throw std::invalid_argument("Sampling grid exceeds 2000000 points");
    const auto weights=weightsOrUnit(q.weights,q.stateDimension);
    SampledContractionReport out;const double low=-std::numeric_limits<double>::infinity();
    out.worst={low,low,low};std::vector<double> worst2,worst1,worstInf;bool invariant=true;
    std::vector<std::size_t> index(D,0);std::vector<double> x(D);
    const auto axisCount=[&](std::size_t i){return q.lower[i]==q.upper[i]?std::size_t(1):q.samplesPerAxis;};
    while(true){
        for(std::size_t i=0;i<D;++i){const auto c=axisCount(i);x[i]=c==1?q.lower[i]:q.lower[i]+(q.upper[i]-q.lower[i])*double(index[i])/double(c-1);}
        const auto J=q.jacobian(x);if(J.rows!=q.stateDimension||J.cols!=q.stateDimension)throw std::invalid_argument("Jacobian size mismatch");
        const auto n=logarithmicNorms(J,weights);++out.samples;
        if(n.two>out.worst.two){out.worst.two=n.two;worst2=x;}
        if(n.one>out.worst.one){out.worst.one=n.one;worst1=x;}
        if(n.infinity>out.worst.infinity){out.worst.infinity=n.infinity;worstInf=x;}
        if(q.field){const auto f=q.field(x);if(f.size()!=q.stateDimension)throw std::invalid_argument("Field size mismatch");
            for(std::size_t i=0;i<q.stateDimension;++i){
                if(axisCount(i)==1)continue;
                if((index[i]==0&&f[i]<0)||(index[i]+1==axisCount(i)&&f[i]>0)){invariant=false;++out.invarianceViolations;}
            }}
        std::size_t a=0;for(;a<D;++a){if(++index[a]<axisCount(a))break;index[a]=0;}
        if(a==D)break;
    }
    const auto [best,name]=bestOf(out.worst);out.rateEstimate=best;out.bestNorm=name;
    out.worstPoint=out.bestNorm=="two"?worst2:out.bestNorm=="one"?worst1:worstInf;
    out.contractiveAtSamples=best<0;
    out.warnings.push_back("SAMPLED: the bound holds at grid points only and is not a certificate over the box");
    if(q.field){out.invariantAtSamples=invariant;
        if(!invariant)out.warnings.push_back("Box is not forward invariant at the samples: contraction inside it does not give a unique steady response");}
    else out.warnings.push_back("Forward invariance of the box was not tested: supply the field, or establish it separately");
    if(!out.contractiveAtSamples)out.warnings.push_back("No negative log norm at the samples in these coordinates: this does not prove the system is non-contracting");
    return out;
}
IntervalContractionReport intervalContraction(const std::vector<double>& lower,const std::vector<double>& upper,std::size_t n,const std::vector<double>& weights){
    if(n<1||n>256||lower.size()!=n*n||upper.size()!=n*n)throw std::invalid_argument("Interval bounds need n*n entries");
    for(std::size_t k=0;k<n*n;++k)if(!std::isfinite(lower[k])||!std::isfinite(upper[k])||lower[k]>upper[k])throw std::invalid_argument("Invalid interval entry");
    const auto w=weightsOrUnit(weights,n);IntervalContractionReport out;ComplexMatrix mid(n,n);double radius=0;
    const auto scaled=[&](const std::vector<double>& m,std::size_t i,std::size_t j){return m[i*n+j]*w[i]/w[j];};
    out.upperBound.one=out.upperBound.infinity=-std::numeric_limits<double>::infinity();
    for(std::size_t i=0;i<n;++i){double row=scaled(upper,i,i),column=scaled(upper,i,i);
        for(std::size_t j=0;j<n;++j){
            const double lo=scaled(lower,i,j),hi=scaled(upper,i,j);mid(i,j)=.5*(lo+hi);radius+=.25*(hi-lo)*(hi-lo);
            if(j!=i){row+=std::max(std::abs(lo),std::abs(hi));column+=std::max(std::abs(scaled(lower,j,i)),std::abs(scaled(upper,j,i)));}
        }
        out.upperBound.infinity=std::max(out.upperBound.infinity,row);out.upperBound.one=std::max(out.upperBound.one,column);}
    out.upperBound.two=largestHermitianEigenvalue(mid)+std::sqrt(radius);
    const auto [best,name]=bestOf(out.upperBound);out.rateBound=best;out.bestNorm=name;out.contractiveForInterval=best<0;
    out.scope="Valid for every Jacobian inside the caller's enclosure, which is trusted and not verified here. No outward rounding. Forward invariance of the set must be established separately";
    return out;
}
}
