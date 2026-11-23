#include <dynsyslab/math/OrbitGramian.hpp>
#include <dynsyslab/math/OrbitIntegrator.hpp>
#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>
#include <string>
namespace dynsyslab::math {
namespace {
using C=std::complex<double>;using V=std::vector<double>;
// (exp(z tau)-1)/z with its removable singularity at z=0.
C phi1(C z,double tau){
    const C x=z*tau;
    if(std::abs(x)<1e-3)return tau*(1.+x/2.+x*x/6.+x*x*x/24.+x*x*x*x/120.);
    return (std::exp(x)-1.)/z;
}
// Cyclic Jacobi for a real symmetric matrix (Golub & Van Loan, Matrix Computations, 8.5).
void jacobi(V& a,V& vectors,std::size_t n){
    vectors.assign(n*n,0);for(std::size_t i=0;i<n;++i)vectors[i*n+i]=1;
    for(int sweep=0;sweep<80;++sweep){
        double off=0,diag=0;
        for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j)(i==j?diag:off)+=a[i*n+j]*a[i*n+j];
        if(off<=1e-30*diag||off==0)return;
        for(std::size_t p=0;p+1<n;++p)for(std::size_t q=p+1;q<n;++q){
            const double apq=a[p*n+q];if(apq==0)continue;
            const double theta=(a[q*n+q]-a[p*n+p])/(2*apq);
            const double t=(theta>=0?1.:-1.)/(std::abs(theta)+std::sqrt(theta*theta+1));
            const double c=1/std::sqrt(t*t+1),s=t*c;
            for(std::size_t k=0;k<n;++k){const double x=a[k*n+p],y=a[k*n+q];a[k*n+p]=c*x-s*y;a[k*n+q]=s*x+c*y;}
            for(std::size_t k=0;k<n;++k){const double x=a[p*n+k],y=a[q*n+k];a[p*n+k]=c*x-s*y;a[q*n+k]=s*x+c*y;}
            for(std::size_t k=0;k<n;++k){const double x=vectors[k*n+p],y=vectors[k*n+q];vectors[k*n+p]=c*x-s*y;vectors[k*n+q]=s*x+c*y;}
        }
    }
}
// Hermitian H=X+iY embeds as the real symmetric [[X,-Y],[Y,X]], doubling each eigenvalue.
V embed(const ComplexVector& h,std::size_t n){
    V a(4*n*n);const auto m=2*n;
    for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j){
        const C z=.5*(h[i*n+j]+std::conj(h[j*n+i]));
        a[i*m+j]=z.real();a[i*m+n+j]=-z.imag();a[(n+i)*m+j]=z.imag();a[(n+i)*m+n+j]=z.real();
    }
    return a;
}
std::vector<double> hankelValues(const ComplexVector& wc,const ComplexVector& wo,std::size_t n){
    const auto m=2*n;V observability=embed(wo,n),vectors;jacobi(observability,vectors,m);
    V root(m*m);
    for(std::size_t i=0;i<m;++i)for(std::size_t j=0;j<m;++j){double sum=0;
        for(std::size_t k=0;k<m;++k)sum+=vectors[i*m+k]*std::sqrt(std::max(observability[k*m+k],0.))*vectors[j*m+k];
        root[i*m+j]=sum;}
    const V control=embed(wc,n);V left(m*m),product(m*m);
    for(std::size_t i=0;i<m;++i)for(std::size_t k=0;k<m;++k){const double r=root[i*m+k];if(r==0)continue;
        for(std::size_t j=0;j<m;++j)left[i*m+j]+=r*control[k*m+j];}
    for(std::size_t i=0;i<m;++i)for(std::size_t k=0;k<m;++k){const double r=left[i*m+k];if(r==0)continue;
        for(std::size_t j=0;j<m;++j)product[i*m+j]+=r*root[k*m+j];}
    for(std::size_t i=0;i<m;++i)for(std::size_t j=i+1;j<m;++j)product[i*m+j]=product[j*m+i]=.5*(product[i*m+j]+product[j*m+i]);
    V unused;jacobi(product,unused,m);V values(m);for(std::size_t i=0;i<m;++i)values[i]=product[i*m+i];
    std::sort(values.rbegin(),values.rend());V out(n);
    for(std::size_t i=0;i<n;++i)out[i]=std::sqrt(std::max(.5*(values[2*i]+values[2*i+1]),0.));
    return out;
}
void validateChannels(const std::vector<ComplexVector>& rows,std::size_t n,const char* name){
    if(rows.empty())throw std::invalid_argument(std::string(name)+" needs at least one channel");
    for(const auto& row:rows){if(row.size()!=n)throw std::invalid_argument(std::string(name)+" channel size differs from the mode count");
        for(C z:row)if(!std::isfinite(z.real())||!std::isfinite(z.imag()))throw std::invalid_argument(std::string(name)+" must be finite");}
}
void requireFinite(double v,const char* m){if(!std::isfinite(v))throw std::invalid_argument(m);}
const JordanMode* findMode(const std::vector<JordanMode>& modes,int n){
    for(const auto& m:modes)if(m.harmonic==n)return &m;return nullptr;
}
void validateModes(const std::vector<JordanMode>& modes){
    std::set<int> seen;
    for(const auto& m:modes){
        if(!seen.insert(m.harmonic).second)throw std::invalid_argument("Duplicate Jordan harmonic");
        for(double v:{m.amplitude.real(),m.amplitude.imag(),m.slope.real(),m.slope.imag()})requireFinite(v,"Jordan coefficients must be finite");
    }
}
}
ModalGramians modalGramians(const ComplexVector& lambda,const std::vector<ComplexVector>& inputs,
    const std::vector<ComplexVector>& outputs,std::optional<double> horizon){
    const auto n=lambda.size();if(n==0||n>128)throw std::invalid_argument("Modal Gramian needs 1..128 modes");
    validateChannels(inputs,n,"Inputs");validateChannels(outputs,n,"Outputs");
    for(C z:lambda){requireFinite(z.real(),"Eigenvalues must be finite");requireFinite(z.imag(),"Eigenvalues must be finite");
        if(!horizon&&!(z.real()<0))throw std::invalid_argument("Infinite-horizon Gramian needs Re(lambda)<0");}
    if(horizon&&(!(*horizon>0)||!std::isfinite(*horizon)))throw std::invalid_argument("Gramian horizon must be positive");
    ModalGramians r;r.dimension=n;r.controllability.resize(n*n);r.observability.resize(n*n);
    const auto kernel=[&](C z){return horizon?phi1(z,*horizon):-1./z;};
    for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j){
        C bb=0,cc=0;for(const auto& b:inputs)bb+=b[i]*std::conj(b[j]);for(const auto& c:outputs)cc+=std::conj(c[i])*c[j];
        r.controllability[i*n+j]=bb*kernel(lambda[i]+std::conj(lambda[j]));
        r.observability[i*n+j]=cc*kernel(std::conj(lambda[i])+lambda[j]);
    }
    double h2=0;
    for(const auto& c:outputs)for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j)
        h2+=(c[i]*r.controllability[i*n+j]*std::conj(c[j])).real();
    r.h2Norm=std::sqrt(std::max(h2,0.));r.hankelSingularValues=hankelValues(r.controllability,r.observability,n);return r;
}
std::array<double,2> oscillatorHankel(double a,double d){
    if(!(a>0)||!(d>0)||!std::isfinite(a)||!std::isfinite(d))throw std::invalid_argument("Oscillator needs a,d>0");
    const double root=std::sqrt(4*a+d*d),denominator=4*a*d;return {(root+d)/denominator,(root-d)/denominator};
}
std::array<double,2> hankelFromSymmetric(const std::array<double,3>& c,const std::array<double,3>& o){
    // Entries of Wc*Wo: trace and determinant give the two squared singular values.
    const double m00=c[0]*o[0]+c[1]*o[1],m01=c[0]*o[1]+c[1]*o[2],m10=c[1]*o[0]+c[2]*o[1],m11=c[1]*o[1]+c[2]*o[2];
    const double half=.5*(m00+m11),gap=.5*(m00-m11),root=std::sqrt(std::max(gap*gap+m01*m10,0.));
    const double upper=half+root,determinant=m00*m11-m01*m10;
    const double lower=upper>0?determinant/upper:0;  // avoids cancellation in half-root
    return {std::sqrt(std::max(upper,0.)),std::sqrt(std::max(lower,0.))};
}
OrbitGramianPair undampedOrbitGramians(const OrbitGramianInput& g){
    for(double v:{g.inputVariance,g.outputVariance,g.frequencySlope,g.horizon,g.inputCross,g.inputPhase,g.outputCross,g.outputEnergy})
        requireFinite(v,"Orbit Gramian inputs must be finite");
    if(!(g.inputVariance>=0)||!(g.outputVariance>=0)||!(g.horizon>0)||!(g.inputPhase>=0)||!(g.outputEnergy>=0))
        throw std::invalid_argument("Orbit Gramian needs nonnegative variances and a positive horizon");
    const double t=g.horizon,s=g.frequencySlope,B=g.inputVariance,C_=g.outputVariance;
    OrbitGramianPair r;
    r.controllability={t*B,s*B*t*t/2+g.inputCross*t,s*s*B*t*t*t/3+s*g.inputCross*t*t+g.inputPhase*t};
    r.observability={s*s*C_*t*t*t/3+s*g.outputCross*t*t+g.outputEnergy*t,s*C_*t*t/2+g.outputCross*t,t*C_};
    r.hankel=hankelFromSymmetric(r.controllability,r.observability);return r;
}
std::array<double,2> collocatedOrbitHankel(double B,double twist,double tau){
    if(!(B>=0)||!(tau>0)||!std::isfinite(B)||!std::isfinite(twist)||!std::isfinite(tau))
        throw std::invalid_argument("Collocated law needs B>=0 and tau>0");
    const double scale=B*std::abs(twist)*tau*tau,root=1/std::sqrt(3.);return {scale*(root+.5),scale*(root-.5)};
}
DampedGramianResult dampedOrbitGramians(const DampedGramianInput& q){
    requireIncreasing(q.energy,"Damped Gramian energy");const auto n=q.energy.size();
    if(q.inputVariance.size()!=n||q.outputVariance.size()!=n||q.damping.size()!=n||q.frequency.size()!=n)
        throw std::invalid_argument("Damped Gramian table sizes differ");
    if(!(q.dampingScale>0)||!(q.initialEnergy>q.targetEnergy)||!(q.targetEnergy>q.finalEnergy))
        throw std::invalid_argument("Damped Gramian needs d>0 and E0>Et>Ef");
    if(q.finalEnergy<q.energy.front()||q.initialEnergy>q.energy.back())throw std::invalid_argument("Energy limits lie outside the table");
    const auto at=[&](const V& f,double e){return interpolateTable(q.energy,f,e);};
    const double Dt=at(q.damping,q.targetEnergy),Ot=at(q.frequency,q.targetEnergy),d=q.dampingScale;
    if(!(Dt>0))throw std::invalid_argument("Averaged damping must be positive at the target");
    static const std::array<double,3> gx{-.7745966692414834,0.,.7745966692414834},gw{5./9.,8./9.,5./9.};
    // Integrates fn(B,C,D,Omega) over [a,b], breaking panels at the table nodes.
    const auto integrate=[&](double a,double b,auto fn){
        double sum=0,left=a;std::size_t i=std::size_t(std::upper_bound(q.energy.begin(),q.energy.end(),a)-q.energy.begin());
        while(left<b){
            const double right=i<n?std::min(q.energy[i],b):b;
            if(right>left){const double mid=.5*(left+right),half=.5*(right-left);
                for(int k=0;k<3;++k){const double e=mid+half*gx[k],D=at(q.damping,e);
                    if(!(D>0))throw std::invalid_argument("Averaged damping must be positive on the quadrature interval");
                    sum+=gw[k]*half*fn(at(q.inputVariance,e),at(q.outputVariance,e),D,at(q.frequency,e));}}
            left=right;++i;
        }
        return sum;};
    DampedGramianResult r;r.targetDamping=Dt;r.targetFrequency=Ot;const double E0=q.initialEnergy,Et=q.targetEnergy,Ef=q.finalEnergy;
    r.cEE=Dt*Dt/d*integrate(Et,E0,[&](double B,double,double D,double){return B/(D*D*D);});
    r.cEphi=Dt/(d*d)*integrate(Et,E0,[&](double B,double,double D,double O){return B*(O-Ot)/(D*D*D);});
    r.cPhiphi=1/(d*d*d)*integrate(Et,E0,[&](double B,double,double D,double O){return B*(O-Ot)*(O-Ot)/(D*D*D);});
    r.oPhiphi=1/d*integrate(Ef,Et,[&](double,double C_,double D,double){return C_/D;});
    r.oEphi=1/(d*d*Dt)*integrate(Ef,Et,[&](double,double C_,double D,double O){return C_*(Ot-O)/D;});
    r.oEE=1/(d*d*d*Dt*Dt)*integrate(Ef,Et,[&](double,double C_,double D,double O){return C_*(Ot-O)*(Ot-O)/D;});
    r.hankel=hankelFromSymmetric({r.cEE,r.cEphi,r.cPhiphi},{r.oEE,r.oEphi,r.oPhiphi});return r;
}
ComplexVector jordanImpulseResponse(const std::vector<JordanMode>& modes,double w,double twist,const V& times){
    validateModes(modes);requireFinite(w,"Frequency must be finite");requireFinite(twist,"Twist must be finite");
    for(const auto& m:modes){const auto* other=findMode(modes,-m.harmonic);
        if(!other||std::abs(other->amplitude-std::conj(m.amplitude))>1e-12*(1+std::abs(m.amplitude))
            ||std::abs(other->slope-std::conj(m.slope))>1e-12*(1+std::abs(m.slope)))
            throw std::invalid_argument("Jordan response needs W_{-n}=conj(W_n)");}
    ComplexVector out;out.reserve(times.size());
    for(double t:times){requireFinite(t,"Response time must be finite");C sum=0;
        for(const auto& m:modes){const double n=m.harmonic,power=std::norm(m.amplitude),rate=2*(std::conj(m.amplitude)*m.slope).real();
            sum+=(C(0,n*rate)-n*n*twist*t*power)*std::exp(C(0,n*w*t));}
        out.push_back(sum);}
    return out;
}
std::map<int,C> frozenKromResidues(const std::vector<JordanMode>& modes){
    validateModes(modes);std::map<int,C> out;
    for(const auto& m:modes){const auto* other=findMode(modes,-m.harmonic);if(other)out[m.harmonic]=C(0,double(m.harmonic))*m.amplitude*other->slope;}
    return out;
}
std::map<int,C> exactInputLift(int n,const std::vector<JordanMode>& modes){
    validateModes(modes);std::map<int,C> out;
    for(const auto& m:modes)out[n+m.harmonic]+=C(0,double(n))*m.slope;
    return out;
}
}
