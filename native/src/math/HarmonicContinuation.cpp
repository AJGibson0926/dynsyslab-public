#include <dynsyslab/math/HarmonicContinuation.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace dynsyslab::math {
namespace {
using C=std::complex<double>;
using V=std::vector<double>;
const double pi=std::acos(-1.);
struct NonFinite {};
double dot(const V& a,const V& b){double s=0;for(std::size_t k=0;k<a.size();++k)s+=a[k]*b[k];return s;}
double maxNorm(const V& a){double s=0;for(double x:a)s=std::max(s,std::abs(x));return s;}
class Balance {
  public:
    explicit Balance(const ContinuationProblem& p):p_(p),d_(p.dimension),K_(std::size_t(p.harmonics)),
        blocks_(2*K_+1),M_(std::max<std::size_t>(64,8*K_+8)),cosine_(M_*(K_+1)),sine_(M_*(K_+1)){
        for(std::size_t j=0;j<M_;++j)for(std::size_t k=0;k<=K_;++k){const double a=2*pi*double(j*k)/double(M_);cosine_[j*(K_+1)+k]=std::cos(a);sine_[j*(K_+1)+k]=std::sin(a);}
    }
    std::size_t unknowns() const {return d_*blocks_;}
    void parameters(double value,double& w,double& a,double& q) const {
        w=p_.frequency;a=p_.amplitude;q=p_.parameter;
        if(p_.variable==ContinuationVariable::Frequency)w=value;else if(p_.variable==ContinuationVariable::Amplitude)a=value;else q=value;
    }
    // Y=[coefficients;value].
    V residual(const V& Y) const {
        double w,a,q;parameters(Y.back(),w,a,q);const auto n=unknowns();V r(n,0.),x(d_),f(d_);
        for(std::size_t j=0;j<M_;++j){
            for(std::size_t s=0;s<d_;++s){const double* c=&Y[s*blocks_];double v=c[0];
                for(std::size_t k=1;k<=K_;++k)v+=c[2*k-1]*cosine_[j*(K_+1)+k]+c[2*k]*sine_[j*(K_+1)+k];x[s]=v;}
            p_.field(x.data(),a*cosine_[j*(K_+1)+1],q,f.data());
            for(std::size_t s=0;s<d_;++s){if(!std::isfinite(f[s]))throw NonFinite{};
                double* R=&r[s*blocks_];R[0]-=f[s]/double(M_);
                for(std::size_t k=1;k<=K_;++k){R[2*k-1]-=2*f[s]*cosine_[j*(K_+1)+k]/double(M_);R[2*k]-=2*f[s]*sine_[j*(K_+1)+k]/double(M_);}}
        }
        for(std::size_t s=0;s<d_;++s)for(std::size_t k=1;k<=K_;++k){r[s*blocks_+2*k-1]+=double(k)*w*Y[s*blocks_+2*k];r[s*blocks_+2*k]-=double(k)*w*Y[s*blocks_+2*k-1];}
        return r;
    }
    // Forward-difference Jacobian n by (n+1).
    std::vector<V> jacobian(const V& Y,const V& r0) const {
        std::vector<V> columns(Y.size());V y=Y;
        for(std::size_t c=0;c<Y.size();++c){const double h=1e-7*(1+std::abs(Y[c]));y[c]=Y[c]+h;auto r=residual(y);y[c]=Y[c];
            for(std::size_t i=0;i<r.size();++i)r[i]=(r[i]-r0[i])/h;columns[c]=std::move(r);}
        return columns;
    }
    std::map<int,C> response(const V& Y) const {
        std::map<int,C> out;const double* c=&Y[p_.outputState*blocks_];out[0]=c[0];
        for(std::size_t k=1;k<=K_;++k)out[int(k)]=.5*C(c[2*k-1],-c[2*k]);
        return out;
    }
    // Spectral radius of the monodromy of the variational equation (RK4, frozen coefficients on stages).
    double floquetRadius(const V& Y) const {
        double w,a,q;parameters(Y.back(),w,a,q);const std::size_t steps=std::max<std::size_t>(256,32*K_);const double T=2*pi/w,h=T/double(steps);
        const auto rate=[&](double t,const std::vector<double>& phi){
            V x(d_),f0(d_),f1(d_);for(std::size_t s=0;s<d_;++s){const double* c=&Y[s*blocks_];double v=c[0];
                for(std::size_t k=1;k<=K_;++k)v+=c[2*k-1]*std::cos(double(k)*w*t)+c[2*k]*std::sin(double(k)*w*t);x[s]=v;}
            const double u=a*std::cos(w*t);p_.field(x.data(),u,q,f0.data());V J(d_*d_);
            for(std::size_t c=0;c<d_;++c){const double step=1e-7*(1+std::abs(x[c])),keep=x[c];x[c]=keep+step;p_.field(x.data(),u,q,f1.data());x[c]=keep;
                for(std::size_t i=0;i<d_;++i){if(!std::isfinite(f1[i]))throw NonFinite{};J[i*d_+c]=(f1[i]-f0[i])/step;}}
            V out(d_*d_,0.);for(std::size_t i=0;i<d_;++i)for(std::size_t j=0;j<d_;++j)for(std::size_t k=0;k<d_;++k)out[i*d_+j]+=J[i*d_+k]*phi[k*d_+j];
            return out;};
        V phi(d_*d_,0.);for(std::size_t i=0;i<d_;++i)phi[i*d_+i]=1;
        const auto add=[](const V& x,const V& g,double f){V o=x;for(std::size_t k=0;k<o.size();++k)o[k]+=f*g[k];return o;};
        for(std::size_t n=0;n<steps;++n){const double t=double(n)*h;
            const auto k1=rate(t,phi),k2=rate(t+.5*h,add(phi,k1,.5*h)),k3=rate(t+.5*h,add(phi,k2,.5*h)),k4=rate(t+h,add(phi,k3,h));
            for(std::size_t e=0;e<phi.size();++e)phi[e]+=h/6*(k1[e]+2*k2[e]+2*k3[e]+k4[e]);}
        ComplexMatrix m(d_,d_);for(std::size_t e=0;e<phi.size();++e)m.data[e]=phi[e];
        double radius=0;for(const auto& z:complexEigenvalues(m))radius=std::max(radius,std::abs(z));
        return radius;
    }
  private:
    const ContinuationProblem& p_;
    std::size_t d_,K_,blocks_,M_;
    V cosine_,sine_;
};
// Solves the bordered system [J;row] x=rhs.
bool solveBordered(const std::vector<V>& columns,const V& row,const V& rhs,V& x){
    const auto n=columns.size();ComplexMatrix m(n,n);
    for(std::size_t c=0;c<n;++c){for(std::size_t i=0;i+1<n;++i)m(i,c)=columns[c][i];m(n-1,c)=row[c];}
    ComplexLu lu(std::move(m));if(lu.singular())return false;
    const auto z=lu.solve(ComplexVector(rhs.begin(),rhs.end()));x.resize(n);
    for(std::size_t k=0;k<n;++k){x[k]=z[k].real();if(!std::isfinite(x[k]))return false;}
    return true;
}
bool tangentAt(const Balance& b,const V& Y,const V& guide,V& tangent){
    const auto J=b.jacobian(Y,b.residual(Y));V rhs(Y.size(),0.);rhs.back()=1;
    if(!solveBordered(J,guide,rhs,tangent))return false;
    const double n=std::sqrt(dot(tangent,tangent));if(!(n>0)||!std::isfinite(n))return false;
    for(auto& x:tangent)x/=n;
    return true;
}
// Newton on [R(Y)=0; row.(Y-anchor)=0].
bool correct(const Balance& b,V& Y,const V& row,const V& anchor,double tolerance,double& residual){
    for(int it=0;it<12;++it){
        const auto r=b.residual(Y);residual=maxNorm(r);
        V rhs(Y.size());for(std::size_t i=0;i<r.size();++i)rhs[i]=-r[i];
        double constraint=0;for(std::size_t k=0;k<Y.size();++k)constraint+=row[k]*(Y[k]-anchor[k]);rhs.back()=-constraint;
        if(residual<=tolerance&&std::abs(constraint)<=1e-9)return true;
        V delta;if(!solveBordered(b.jacobian(Y,r),row,rhs,delta))return false;
        for(std::size_t k=0;k<Y.size();++k)Y[k]+=delta[k];
    }
    residual=maxNorm(b.residual(Y));
    return residual<=tolerance;
}
void validate(const ContinuationProblem& p){
    if(p.dimension<1||p.dimension>32||!p.field)throw std::invalid_argument("Continuation needs a field of dimension 1..32");
    if(p.harmonics<1||p.harmonics>64)throw std::invalid_argument("Harmonic count must be 1..64");
    if(p.dimension*(2*std::size_t(p.harmonics)+1)>1024)throw std::invalid_argument("Harmonic balance system exceeds 1024 unknowns");
    for(double x:{p.frequency,p.amplitude,p.parameter,p.lower,p.upper})if(!std::isfinite(x))throw std::invalid_argument("Continuation values must be finite");
    if(!(p.lower<p.upper))throw std::invalid_argument("Continuation needs lower<upper");
    if(!(p.minimumStep>0)||!(p.step>=p.minimumStep)||!(p.maximumStep>=p.step)||!std::isfinite(p.maximumStep))throw std::invalid_argument("Need 0<minimum_step<=step<=maximum_step");
    if(p.maximumSteps<1||p.maximumSteps>100000||!(p.tolerance>0)||(p.direction!=1&&p.direction!=-1))throw std::invalid_argument("Invalid continuation controls");
    if(p.outputState>=p.dimension)throw std::invalid_argument("Output state out of range");
    if(!p.initialCoefficients.empty()&&p.initialCoefficients.size()!=p.dimension*(2*std::size_t(p.harmonics)+1))throw std::invalid_argument("Initial coefficients need d(2K+1) entries");
    for(double x:p.initialCoefficients)if(!std::isfinite(x))throw std::invalid_argument("Initial coefficients must be finite");
    const double start=p.variable==ContinuationVariable::Frequency?p.frequency:p.variable==ContinuationVariable::Amplitude?p.amplitude:p.parameter;
    if(start<p.lower||start>p.upper)throw std::invalid_argument("Start value lies outside the continuation bounds");
    if(p.variable!=ContinuationVariable::Frequency&&!(p.frequency>0))throw std::invalid_argument("Forcing frequency must be positive");
    if(p.variable==ContinuationVariable::Frequency&&!(p.lower>0))throw std::invalid_argument("Frequency continuation needs a positive lower bound");
}
}
ContinuationResult continueHarmonicBalance(const ContinuationProblem& p){
    validate(p);const Balance b(p);const auto n=b.unknowns();ContinuationResult out;
    V Y(n+1,0.);if(!p.initialCoefficients.empty())std::copy(p.initialCoefficients.begin(),p.initialCoefficients.end(),Y.begin());
    Y[n]=p.variable==ContinuationVariable::Frequency?p.frequency:p.variable==ContinuationVariable::Amplitude?p.amplitude:p.parameter;
    const auto record=[&](const V& y,double residual){
        ContinuationPoint q;q.value=y[n];q.residual=residual;q.coefficients.assign(y.begin(),y.begin()+std::ptrdiff_t(n));q.response=b.response(y);
        if(p.computeStability){const double r=b.floquetRadius(y);q.floquetRadius=r;q.stable=r<1;}
        out.points.push_back(std::move(q));};
    const auto stop=[&](const char* reason,bool lost,const std::string& note){out.termination=reason;out.branchLost=lost;out.warnings.push_back(note);};
    try{
        V fix(n+1,0.);fix[n]=1;double residual=0;const V anchor=Y;
        if(!correct(b,Y,fix,anchor,p.tolerance,residual)){stop("corrector_failed",true,"No periodic response found at the start value");return out;}
        record(Y,residual);
        V tangent;fix[n]=double(p.direction);
        if(!tangentAt(b,Y,fix,tangent)){stop("singular_tangent",true,"Singular Jacobian at the start: the branch cannot be followed from here");return out;}
        if(tangent[n]*p.direction<0)for(auto& x:tangent)x=-x;
        double ds=p.step;
        for(std::size_t step=0;step<p.maximumSteps;++step){
            V next=Y;for(std::size_t k=0;k<=n;++k)next[k]+=ds*tangent[k];const V predicted=next;
            V newTangent;bool ok=correct(b,next,tangent,predicted,p.tolerance,residual)&&tangentAt(b,next,tangent,newTangent);
            const double turning=ok?dot(tangent,newTangent):-1;
            if(!ok||turning<.8){
                ds*=.5;
                if(ds<p.minimumStep){
                    if(ok)stop("tangent_jump",true,"Tangent turned abruptly at the minimum step: possible branch point or a jump to another branch");
                    else stop("corrector_failed",true,"Corrector failed at the minimum step: the branch is lost here");
                    return out;
                }
                continue;
            }
            if(next[n]<p.lower||next[n]>p.upper){stop("parameter_bound",false,"Reached the continuation bound");return out;}
            if(p.variable==ContinuationVariable::Frequency&&!(next[n]>0)){stop("parameter_bound",false,"Frequency reached zero");return out;}
            if(tangent[n]*newTangent[n]<0){
                const double theta=tangent[n]/(tangent[n]-newTangent[n]);
                const auto before=b.response(Y).at(1),after=b.response(next).at(1);
                out.folds.push_back({Y[n]+ds*tangent[n]*theta/2,out.points.size()-1,(1-theta)*std::abs(before)+theta*std::abs(after)});
            }
            Y=std::move(next);tangent=std::move(newTangent);record(Y,residual);
            if(turning>.98)ds=std::min(ds*1.3,p.maximumStep);
        }
        stop("step_limit",false,"Step limit reached before the continuation bound");
    }catch(const NonFinite&){stop("non_finite_field",true,"The field returned a non-finite value: the branch is lost here");}
    return out;
}
}
