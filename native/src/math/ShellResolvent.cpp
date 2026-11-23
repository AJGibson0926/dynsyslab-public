#include <dynsyslab/math/DampedOrbit.hpp>
#include <dynsyslab/math/ShellResolvent.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <set>
#include <stdexcept>
namespace dynsyslab::math {
namespace {
using C=std::complex<double>;
using V=std::vector<double>;
const double pi=std::acos(-1.);
// Eight-point Gauss-Legendre rule on [-1,1] (Abramowitz & Stegun, Table 25.4).
const std::array<double,4> gaussX{.1834346424956498,.5255324099163290,.7966664774136267,.9602898564975363};
const std::array<double,4> gaussW{.3626837833783620,.3137066458778873,.2223810344533745,.1012285362903763};
V realDerivative(const V& x,const V& y){
    ComplexVector z(y.begin(),y.end());const auto d=hermiteNodeDerivatives(x,z);
    V out(d.size());for(std::size_t k=0;k<d.size();++k)out[k]=d[k].real();return out;
}
struct Hermite {
    const V& x;const V& y;V slope;
    Hermite(const V& nodes,const V& values):x(nodes),y(values),slope(realDerivative(nodes,values)){}
    std::size_t panel(double e) const {
        auto it=std::upper_bound(x.begin(),x.end(),e);
        if(it==x.begin())return 0;
        const auto k=std::size_t(it-x.begin())-1;return std::min(k,x.size()-2);
    }
    double value(double e) const {const auto k=panel(e);const double h=x[k+1]-x[k],u=(e-x[k])/h;
        const double h00=(1+2*u)*(1-u)*(1-u),h10=u*(1-u)*(1-u),h01=u*u*(3-2*u),h11=u*u*(u-1);
        return h00*y[k]+h10*h*slope[k]+h01*y[k+1]+h11*h*slope[k+1];}
    double derivative(double e) const {const auto k=panel(e);const double h=x[k+1]-x[k],u=(e-x[k])/h;
        const double d00=6*u*(u-1)/h,d10=(1-u)*(1-3*u),d01=-d00,d11=u*(3*u-2);
        return d00*y[k]+d10*slope[k]+d01*y[k+1]+d11*slope[k+1];}
};
V harmonicWeight(const ShellResolventTable& t,const ComplexVector& g){
    V q(g.size());for(std::size_t k=0;k<q.size();++k)q[k]=t.density[k]*std::norm(g[k]);return q;
}
// int_0^L (qa+qs u)/(da+beta u) du, with da+beta u off the origin on the panel.
C panelRationalIntegral(C qa,C qs,C da,C beta,double L){
    const C x=beta*L/da;
    if(std::abs(x)<1e-3){
        const C a1=1.-x/2.+x*x/3.-x*x*x/4.,a2=.5-x/3.+x*x/4.-x*x*x/5.;
        return (L/da)*(qa*a1+qs*L*a2);
    }
    return qs*L/beta+(qa-qs*da/beta)*(std::log(da+beta*L)-std::log(da))/beta;
}
}
void validateShellTable(const ShellResolventTable& t){
    const auto n=t.energy.size();
    if(n<3||t.frequency.size()!=n||t.density.size()!=n)throw std::invalid_argument("Shell table needs at least three matching nodes");
    if(t.harmonics.empty())throw std::invalid_argument("Shell table needs harmonics");
    for(std::size_t k=0;k<n;++k){
        if(!std::isfinite(t.energy[k])||(k&&t.energy[k]<=t.energy[k-1]))throw std::invalid_argument("Energies must increase strictly");
        if(!std::isfinite(t.frequency[k])||t.frequency[k]<=0)throw std::invalid_argument("Orbit frequencies must be positive");
        if(!std::isfinite(t.density[k])||t.density[k]<0)throw std::invalid_argument("Densities must be nonnegative");
    }
    for(const auto& [h,g]:t.harmonics){(void)h;if(g.size()!=n)throw std::invalid_argument("Harmonic table size mismatch");
        for(const auto& z:g)if(!std::isfinite(z.real())||!std::isfinite(z.imag()))throw std::invalid_argument("Harmonics must be finite");}
}
ShellResolventTable shellTableFromOrbitShells(const std::vector<OrbitShell>& shells){
    if(shells.size()<3)throw std::invalid_argument("Shell table needs at least three shells");
    ShellResolventTable t;const auto n=shells.size();std::set<int> keys;
    for(const auto& s:shells){t.energy.push_back(s.energy);t.frequency.push_back(s.series.angularFrequency);
        for(const auto& [h,c]:s.series.coefficients){(void)c;keys.insert(h);}}
    for(std::size_t k=1;k<n;++k)if(!(t.energy[k]>t.energy[k-1]))throw std::invalid_argument("Shell energies must increase strictly");
    t.density.resize(n);
    for(std::size_t k=0;k<n;++k){
        if(!std::isfinite(shells[k].weight)||shells[k].weight<0)throw std::invalid_argument("Shell weights must be nonnegative");
        const double left=k?t.energy[k]-t.energy[k-1]:0,right=k+1<n?t.energy[k+1]-t.energy[k]:0;
        t.density[k]=shells[k].weight/(.5*(left+right));
    }
    for(int h:keys){ComplexVector g(n);
        for(std::size_t k=0;k<n;++k){const auto it=shells[k].series.coefficients.find(h);if(it!=shells[k].series.coefficients.end())g[k]=it->second;}
        t.harmonics.emplace(h,std::move(g));}
    validateShellTable(t);
    return t;
}
PoleResidueModel orbitHarmonicModel(const FourierSeries& series,double phase){
    if(series.coefficients.empty())throw std::invalid_argument("Orbit series needs coefficients");
    if(!std::isfinite(series.angularFrequency)||!std::isfinite(phase))throw std::invalid_argument("Orbit frequency and phase must be finite");
    PoleResidueModel m;
    for(const auto& [n,g]:series.coefficients)m.poles.push_back({C(0,n*series.angularFrequency),{g*std::exp(C(0,n*phase))}});
    validatePoleModel(m);
    return m;
}
ShellBoundaryValue shellBoundaryValue(const ShellResolventTable& t,double w,const ShellBoundaryOptions& o){
    validateShellTable(t);if(!std::isfinite(w))throw std::invalid_argument("Frequency must be finite");
    if(!(o.slopeTolerance>0))throw std::invalid_argument("Slope tolerance must be positive");
    if(!o.damping.empty()&&(o.damping.size()!=t.energy.size()||!(o.dampingScale>0)))throw std::invalid_argument("Damping table needs matching size and positive scale");
    ShellBoundaryValue out;out.frequency=w;
    const Hermite Om(t.energy,t.frequency);
    const double lo=t.energy.front(),hi=t.energy.back(),edge=1e-12*(hi-lo);
    double principal=0;
    for(const auto& [n,g]:t.harmonics){
        const V qv=harmonicWeight(t,g);const Hermite q(t.energy,qv);
        if(n==0){
            double mass=0;
            for(std::size_t k=0;k+1<t.energy.size();++k){const double c=.5*(t.energy[k]+t.energy[k+1]),r=.5*(t.energy[k+1]-t.energy[k]);
                for(std::size_t a=0;a<4;++a)mass+=r*gaussW[a]*(q.value(c-r*gaussX[a])+q.value(c+r*gaussX[a]));}
            out.zeroHarmonicMass+=mass;
            if(w!=0)principal+=mass/w;
            continue;
        }
        const auto f=[&](double e){return w-n*Om.value(e);};
        struct Root {double e,slope,A;};std::vector<Root> roots;
        for(std::size_t k=0;k+1<t.energy.size();++k){
            double a=t.energy[k],b=t.energy[k+1],fa=f(a),fb=f(b);
            if(fa==0&&k>0)continue; // A node root belongs to the panel on its left.
            if(fa*fb>0)continue;
            for(int it=0;it<200&&(b-a)>1e-15*(1+std::abs(a)+std::abs(b));++it){const double m=.5*(a+b),fm=f(m);
                if(fm==0){a=b=m;break;}
                if(fa*fm<0){b=m;fb=fm;}else{a=m;fa=fm;}}
            double e=.5*(a+b);
            // Snap to a table node so the panel cuts cannot create a sliver around the root.
            if(std::abs(e-t.energy[k])<=edge)e=t.energy[k];else if(std::abs(e-t.energy[k+1])<=edge)e=t.energy[k+1];
            const double slope=Om.derivative(e),ns=std::abs(n*slope);
            if(e<=lo+edge||e>=hi-edge){out.valid=false;out.warnings.push_back("Harmonic "+std::to_string(n)+": root at a table edge is excluded");continue;}
            if(ns<=o.slopeTolerance){out.valid=false;out.warnings.push_back("Harmonic "+std::to_string(n)+": degenerate root (vanishing twist): weight 1/|n Omega'| is not finite");continue;}
            const double qr=std::max(0.,q.value(e));
            ShellBoundaryRoot root;root.harmonic=n;root.energy=e;root.frequencySlope=slope;root.density=qr;root.weight=qr/ns;
            if(!o.damping.empty()){const Hermite D(t.energy,o.damping);root.stationaryPhaseWeight=shellWeight(o.dampingScale,D.value(e),n,slope);}
            out.roots.push_back(root);out.absorptive+=pi*root.weight;
            roots.push_back({e,slope,qr/(-n*slope)});
        }
        // PV int q/(w-n Omega) = int [q/(w-n Omega) - sum A_r/(E-E_r)] + sum A_r log|(hi-E_r)/(E_r-lo)|.
        V cuts=t.energy;for(const auto& r:roots)cuts.push_back(r.e);std::sort(cuts.begin(),cuts.end());
        cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
        double integral=0;
        for(std::size_t k=0;k+1<cuts.size();++k){const double c=.5*(cuts[k]+cuts[k+1]),r=.5*(cuts[k+1]-cuts[k]);
            const auto eval=[&](double e){double s=q.value(e)/f(e);for(const auto& rt:roots)s-=rt.A/(e-rt.e);return s;};
            for(std::size_t a=0;a<4;++a)integral+=r*gaussW[a]*(eval(c-r*gaussX[a])+eval(c+r*gaussX[a]));}
        for(const auto& r:roots)integral+=r.A*std::log((hi-r.e)/(r.e-lo));
        principal+=integral;
    }
    out.dispersive=-principal;out.value=C(out.absorptive,out.dispersive);
    return out;
}
PoleResidueModel jordanOrbitModel(const std::vector<JordanMode>& modes,double w,double twist){
    // Reuses the OrbitGramian validation: real observable, W_{-n}=conj(W_n).
    (void)jordanImpulseResponse(modes,w,twist,{0.});
    PoleResidueModel m;
    for(const auto& mode:modes){const double n=mode.harmonic,power=std::norm(mode.amplitude),rate=2*(std::conj(mode.amplitude)*mode.slope).real();
        m.poles.push_back({C(0,n*w),{C(0,n*rate),C(-n*n*twist*power)}});}
    validatePoleModel(m);
    return m;
}
C dampedShellPairing(const ShellResolventTable& t,C s,const V& rate,double d){
    validateShellTable(t);
    if(rate.size()!=t.energy.size())throw std::invalid_argument("Damping rate table size mismatch");
    if(!(d>0)||!std::isfinite(d)||!(s.real()>=0)||!std::isfinite(s.imag()))throw std::invalid_argument("Damped pairing needs d>0 and Re s>=0");
    for(double g:rate)if(!std::isfinite(g)||g<=0)throw std::invalid_argument("Damping rates must be positive");
    C total{};
    for(const auto& [n,g]:t.harmonics){const V q=harmonicWeight(t,g);
        for(std::size_t k=0;k+1<t.energy.size();++k){
            const double L=t.energy[k+1]-t.energy[k];
            const C da=s-C(0,n*t.frequency[k])+d*rate[k],db=s-C(0,n*t.frequency[k+1])+d*rate[k+1];
            total+=panelRationalIntegral(q[k],(q[k+1]-q[k])/L,da,(db-da)/L,L);
        }
    }
    return total;
}
RationalResolventCombination rationalKernelResolvent(const RationalCircleKernel& kernel,double theta){
    if(!std::isfinite(theta))throw std::invalid_argument("Spectral angle must be finite");
    RationalResolventCombination out;const double eps=kernel.bandwidth();
    for(std::size_t j=0;j<kernel.poles().size();++j){
        const C p=kernel.poles()[j],a=kernel.residues()[j];
        const C z=std::exp(C(0,theta)-C(0,eps)*p),y=std::exp(C(0,theta)-C(0,eps)*std::conj(p));
        // (z+mu)/(z-mu) = -1 + 2 z/(z-mu).
        out.constant+=(-a+std::conj(a))/(4*pi);
        out.terms.push_back({z,a/(2*pi)});out.terms.push_back({y,-std::conj(a)/(2*pi)});
    }
    return out;
}
C applyRationalCombination(const RationalResolventCombination& r,C mu){
    C sum=r.constant;
    for(const auto& t:r.terms){const C d=t.node-mu;if(d==C{})throw std::domain_error("Kernel evaluated at a resolvent node");
        sum+=t.coefficient*t.node/d;}
    return sum;
}
PoleResidueModel rationalProjectedModel(const PoleResidueModel& model,double h,double theta,const RationalCircleKernel& kernel){
    validatePoleModel(model);if(!(h>0)||!std::isfinite(h))throw std::invalid_argument("Positive sample interval required");
    const auto r=rationalKernelResolvent(kernel,theta);PoleResidueModel out;
    for(const auto& p:model.poles){
        if(p.laurent.size()>2)throw std::invalid_argument("Rational projection supports poles of order at most 2");
        const C mu=std::exp(p.eigenvalue*h),value=applyRationalCombination(r,mu);
        ResolventPole q{p.eigenvalue,{value*p.laurent[0]}};
        if(p.laurent.size()==2){
            C slope{};for(const auto& t:r.terms){const C d=t.node-mu;slope+=t.coefficient*t.node/(d*d);}
            q.laurent[0]+=slope*mu*h*p.laurent[1];
            q.laurent.push_back(value*p.laurent[1]);
        }
        out.poles.push_back(std::move(q));
    }
    return out;
}
}
