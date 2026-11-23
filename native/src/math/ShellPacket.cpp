#include <dynsyslab/math/ShellPacket.hpp>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <set>
#include <stdexcept>
// Rational kernels and spectral packets: Colbrook, Drysdale & Horning, Rigged
// DMD, arXiv:2405.00782. Unitary spectral measures on the circle: Colbrook,
// mpEDMD, SIAM J. Numer. Anal. 61 (2023). Pearson correlation is standard.
namespace dynsyslab::math {
namespace {
using C=std::complex<double>;
const double pi=std::acos(-1.);
void positive(double x,const char* message){if(!std::isfinite(x)||x<=0)throw std::invalid_argument(message);}
void finite(double x,const char* message){if(!std::isfinite(x))throw std::invalid_argument(message);}
std::optional<double> pearson(const std::vector<double>& a,const std::vector<double>& b,const std::vector<double>& w){
    double mass=0,ma=0,mb=0;
    for(std::size_t i=0;i<a.size();++i){mass+=w[i];ma+=w[i]*a[i];mb+=w[i]*b[i];}
    if(mass<=0)return std::nullopt;
    ma/=mass;mb/=mass;double aa=0,bb=0,ab=0;
    for(std::size_t i=0;i<a.size();++i){const double x=a[i]-ma,y=b[i]-mb;aa+=w[i]*x*x;bb+=w[i]*y*y;ab+=w[i]*x*y;}
    const double scale=std::max(std::abs(ma),std::abs(mb));
    const double floor=mass*1e-28*scale*scale;
    if(aa<=floor||bb<=floor)return std::nullopt;
    return std::clamp(ab/std::sqrt(aa*bb),-1.,1.);
}
void validateShells(const std::vector<OrbitShell>& shells){
    if(shells.empty()||shells.size()>1000000)throw std::invalid_argument("Shell table needs 1..1000000 shells");
    for(std::size_t i=0;i<shells.size();++i){
        finite(shells[i].energy,"Shell energy must be finite");
        positive(shells[i].series.angularFrequency,"Shell Omega must be positive and finite");
        finite(shells[i].weight,"Shell weight must be finite");
        if(shells[i].weight<0)throw std::invalid_argument("Shell weight must be nonnegative");
        if(i&&shells[i].energy<=shells[i-1].energy)throw std::invalid_argument("Shell energies must increase strictly");
    }
}
}
double kernelHalfWidth(const RationalCircleKernel& kernel){
    const double peak=kernel(0);
    if(!(peak>0))throw std::domain_error("Kernel peak must be positive");
    const double half=peak/2,step=std::min(kernel.bandwidth()/64,pi/1024);
    double lo=0;
    for(double x=step;x<=pi;x+=step){
        if(kernel(x)<=half){
            double hi=x;
            for(int k=0;k<80&&hi-lo>1e-14*(hi+lo);++k){const double m=(lo+hi)/2;(kernel(m)>half?lo:hi)=m;}
            return (lo+hi)/2;
        }
        lo=x;
    }
    return pi;
}
double shellPacketWeight(double alpha,double theta,const RationalCircleKernel& kernel){
    finite(alpha,"Finite angle required");finite(theta,"Finite angle required");
    return kernel(theta-alpha);
}
FourierSeries predictedShellPacket(const FourierSeries& orbit,double dt,double theta,
    const RationalCircleKernel& kernel,double delay){
    finite(delay,"Finite target delay required");
    auto out=spectralPacket(orbit,dt,theta,kernel);
    if(delay!=0)for(auto& [n,c]:out.coefficients)
        c*=std::exp(C(0,std::remainder(n*orbit.angularFrequency*delay*dt,2*pi)));
    return out;
}
ShellPacketPrediction predictShellPackets(const std::vector<OrbitShell>& shells,double dt,double theta,
    const RationalCircleKernel& kernel,double delay){
    validateShells(shells);
    ShellPacketPrediction out;out.halfWidth=kernelHalfWidth(kernel);
    std::vector<WeightedOrbitSpectrum> weighted;weighted.reserve(shells.size());
    for(const auto& shell:shells){
        out.packets.push_back(predictedShellPacket(shell.series,dt,theta,kernel,delay));
        out.shellEnergy.push_back(meanSquare(out.packets.back()));
        weighted.push_back({shell.series,shell.weight});
    }
    const auto density=smoothOrbitSpectra(weighted,dt,{theta},kernel);
    out.density=density.density[0];out.totalMass=density.totalMass;out.zeroHarmonicMass=density.zeroHarmonicMass;
    return out;
}
std::vector<ShellAlias> shellAliasCatalogue(const std::vector<OrbitShell>& shells,double dt,double theta,
    int maximumHarmonic,double halfWidth){
    validateShells(shells);positive(dt,"Positive finite sample interval required");
    finite(theta,"Finite angle required");positive(halfWidth,"Positive finite half-width required");
    if(maximumHarmonic<1||maximumHarmonic>100000)throw std::invalid_argument("Maximum harmonic must be 1..100000");
    if(shells.size()<2)throw std::invalid_argument("Alias catalogue needs at least two shells");
    std::set<int> harmonics;
    for(const auto& s:shells)for(const auto& [n,c]:s.series.coefficients)
        if(n&&std::abs(n)<=maximumHarmonic&&c!=C{})harmonics.insert(n);
    std::vector<ShellAlias> out;
    for(int n:harmonics)for(std::size_t i=0;i+1<shells.size();++i){
        const double e0=shells[i].energy,e1=shells[i+1].energy;
        const double w0=shells[i].series.angularFrequency,w1=shells[i+1].series.angularFrequency;
        const double a=n*w0*dt,b=n*w1*dt,slope=(w1-w0)/(e1-e0);
        const bool last=i+2==shells.size();
        if(a==b){ // Flat segment: a resonance here has no finite ridge width.
            const double k=std::round((a-theta)/(2*pi));
            if(std::abs(a-theta-2*pi*k)<=1e-14*(1+std::abs(a)))
                out.push_back({n,(long long)k,e0,w0,slope,std::nullopt});
            continue;
        }
        const double first=std::ceil((std::min(a,b)-theta)/(2*pi)),final=std::floor((std::max(a,b)-theta)/(2*pi));
        for(double k=first;k<=final;++k){
            const double t=(theta+2*pi*k-a)/(b-a);
            if(t<0||t>1||(t==1&&!last))continue; // Each table node belongs to one segment.
            out.push_back({n,(long long)k,e0+t*(e1-e0),w0+t*(w1-w0),slope,
                std::optional<double>(halfWidth/(std::abs(n)*dt*std::abs(slope)))});
            if(out.size()>2000000)throw std::invalid_argument("Alias catalogue exceeds 2000000 entries");
        }
    }
    return out;
}
ShellPacketScore scoreShellPacket(const std::vector<C>& computed,const std::vector<C>& predicted,
    const std::vector<std::size_t>& shellIndex,const std::vector<double>& weights){
    const auto n=computed.size();
    if(!n||predicted.size()!=n)throw std::invalid_argument("Computed and predicted packets need equal nonzero length");
    if(!shellIndex.empty()&&shellIndex.size()!=n)throw std::invalid_argument("One shell index required per sample");
    if(!weights.empty()&&weights.size()!=n)throw std::invalid_argument("One weight required per sample");
    std::vector<double> w=weights.empty()?std::vector<double>(n,1.):weights;
    for(double x:w)if(!std::isfinite(x)||x<0)throw std::invalid_argument("Weights must be finite and nonnegative");
    double cc=0,pp=0,error=0;C cross{};
    std::vector<double> mc(n),mp(n);
    for(std::size_t i=0;i<n;++i){
        if(!std::isfinite(std::abs(computed[i]))||!std::isfinite(std::abs(predicted[i])))
            throw std::invalid_argument("Packet samples must be finite");
        cc+=w[i]*std::norm(computed[i]);pp+=w[i]*std::norm(predicted[i]);
        error+=w[i]*std::norm(computed[i]-predicted[i]);cross+=w[i]*std::conj(predicted[i])*computed[i];
        mc[i]=std::abs(computed[i]);mp[i]=std::abs(predicted[i]);
    }
    if(!(pp>0)||!(cc>0))throw std::invalid_argument("Packets must have positive weighted norm");
    ShellPacketScore out;
    out.relativeL2Error=std::sqrt(error/pp);out.complexCorrelation=cross/std::sqrt(cc*pp);
    out.complexCorrelationMagnitude=std::abs(out.complexCorrelation);
    out.magnitudeCorrelation=pearson(mc,mp,w);
    std::size_t groups=n;
    if(!shellIndex.empty()){
        groups=*std::max_element(shellIndex.begin(),shellIndex.end())+1;
        if(groups>n)throw std::invalid_argument("Shell indices must be compact");
    }
    std::vector<double> ec(groups),ep(groups),mass(groups);
    for(std::size_t i=0;i<n;++i){
        const auto g=shellIndex.empty()?i:shellIndex[i];
        ec[g]+=w[i]*std::norm(computed[i]);ep[g]+=w[i]*std::norm(predicted[i]);mass[g]=1;
    }
    out.energyProfileCorrelation=pearson(ec,ep,mass);
    return out;
}
DelayTargetDiagnostic diagnoseDelayTarget(long long first,long long last,long long target,
    const RationalCircleKernel& kernel,double tolerance){
    if(first>last)throw std::invalid_argument("Delay window needs first_lag<=last_lag");
    if(!std::isfinite(tolerance)||tolerance<0||tolerance>=1)
        throw std::invalid_argument("Asymmetry tolerance must lie in [0,1)");
    DelayTargetDiagnostic out;
    out.lagsBefore=target-first;out.lagsAfter=last-target;
    const bool inside=out.lagsBefore>=0&&out.lagsAfter>=0;
    out.retainedSymmetricLag=inside?std::min(out.lagsBefore,out.lagsAfter):0;
    const long long span=out.lagsBefore+out.lagsAfter;
    out.asymmetry=inside&&span>0?double(std::llabs(out.lagsBefore-out.lagsAfter))/double(span):1;
    out.symmetricTailBound=kernel.symmetricTimeTailBound(std::size_t(out.retainedSymmetricLag));
    out.oneSided=out.asymmetry>tolerance;
    if(out.oneSided)out.warning="Target is one-sided in its delay window; the symmetric kernel is truncated at lag "+
        std::to_string(out.retainedSymmetricLag)+". Centre the target or treat the packet as biased.";
    return out;
}
ResolvableShellEstimate resolvableShellEstimate(double constant,double scale,double halfWidth){
    positive(constant,"Positive empirical constant required");positive(scale,"Positive angular scale required");
    positive(halfWidth,"Positive half-width required");
    return {constant*scale/halfWidth,true,
        "Empirical heuristic n_last=c*w/h; c is caller-calibrated, w is a caller angular scale; not a resolution guarantee"};
}
}
