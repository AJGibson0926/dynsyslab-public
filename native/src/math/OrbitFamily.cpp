#include <dynsyslab/math/OrbitFamily.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace dynsyslab::math {
namespace {
const double kPi=std::acos(-1.);
using V=std::vector<double>;
double dot(const PlanarPoint& a,const PlanarPoint& b){return a[0]*b[0]+a[1]*b[1];}
double cross(const PlanarPoint& a,const PlanarPoint& b){return a[0]*b[1]-a[1]*b[0];}
PlanarPoint minus(const PlanarPoint& a,const PlanarPoint& b){return {a[0]-b[0],a[1]-b[1]};}
PlanarPoint unit(const PlanarPoint& d){const double n=std::hypot(d[0],d[1]);
    if(!(n>0)||!std::isfinite(n))throw std::invalid_argument("Section direction must be nonzero");return {d[0]/n,d[1]/n};}
void validate(const PlanarOrbitSystem& s){
    if(!s.drift||!s.hamiltonian)throw std::invalid_argument("Orbit system needs drift and Hamiltonian");
    if(!(s.maximumRadius>0)||!std::isfinite(s.maximumRadius))throw std::invalid_argument("Orbit system needs a positive maximum radius");
    if(!s.density&&(s.coordinate<0||s.coordinate>1||s.momentum<0||s.momentum>1||s.coordinate==s.momentum))
        throw std::invalid_argument("Canonical pair must be the two distinct state indices");
}
// Five-point Gauss-Legendre on [0,1] for the Poincare-lemma primitive.
double radialWeight(const PlanarScalar& rho,const PlanarPoint& c,const PlanarPoint& u){
    static const std::array<double,5> x{.04691007703066800,.23076534494715845,.5,.76923465505284155,.95308992296933200};
    static const std::array<double,5> w{.11846344252809454,.23931433524968324,.28444444444444444,.23931433524968324,.11846344252809454};
    double sum=0;for(int i=0;i<5;++i)sum+=w[i]*x[i]*rho({c[0]+x[i]*u[0],c[1]+x[i]*u[1]});return sum;
}
FlowField augmented(const PlanarOrbitSystem& s){
    return [&s](const V& x,V& dx){
        const PlanarPoint p{x[0],x[1]},f=s.drift(p);dx[0]=f[0];dx[1]=f[1];
        if(s.density){const auto u=minus(p,s.centre);dx[2]=cross(u,f)*radialWeight(s.density,s.centre,u);}
        else dx[2]=x[s.momentum]*f[s.coordinate];
    };
}
FlowField planar(const PlanarOrbitSystem& s,double sign){
    return [&s,sign](const V& x,V& dx){const auto f=s.drift({x[0],x[1]});dx[0]=sign*f[0];dx[1]=sign*f[1];};
}
}
PlanarPoint orbitStart(const PlanarOrbitSystem& s,double energy){
    validate(s);if(!std::isfinite(energy))throw std::invalid_argument("Energy must be finite");
    const auto d=unit(s.sectionDirection);
    const auto at=[&](double r){return s.hamiltonian({s.centre[0]+r*d[0],s.centre[1]+r*d[1]})-energy;};
    constexpr int scan=256;double lo=0,flo=at(0),hi=0,fhi=flo;bool found=false;
    if(flo==0)throw std::invalid_argument("Energy equals the centre value");
    for(int j=1;j<=scan&&!found;++j){
        hi=s.maximumRadius*j/scan;fhi=at(hi);
        if(fhi==0||(flo<0)!=(fhi<0))found=true;else{lo=hi;flo=fhi;}
    }
    if(!found)throw std::invalid_argument("Energy level not bracketed on the section ray");
    for(int k=0;k<200&&hi-lo>4*std::numeric_limits<double>::epsilon()*hi;++k){
        const double mid=.5*(lo+hi),fm=at(mid);
        if(fm==0){lo=hi=mid;break;}
        if((fm<0)==(flo<0)){lo=mid;flo=fm;}else{hi=mid;fhi=fm;}
    }
    const double r=.5*(lo+hi);return {s.centre[0]+r*d[0],s.centre[1]+r*d[1]};
}
PeriodicOrbit periodicOrbit(const PlanarOrbitSystem& s,double energy,const OrbitFamilyOptions& o){
    const auto start=orbitStart(s,energy);const auto d=unit(s.sectionDirection);const PlanarPoint n{-d[1],d[0]};
    const auto field=augmented(s);const auto f0=s.drift(start);const double side=dot(n,f0);
    if(side==0||!std::isfinite(side))throw std::invalid_argument("Flow is tangent to the section at the orbit start");
    const double orient=side>0?1:-1;const auto g=[&](const V& x){return orient*dot(n,minus({x[0],x[1]},s.centre));};
    const double radius=std::hypot(start[0]-s.centre[0],start[1]-s.centre[1]);
    V x{start[0],start[1],0},previous,trial;double t=0,h=o.flow.initialStep>0?o.flow.initialStep:1e-2,winding=0;
    for(std::size_t step=0;step<o.flow.maximumSteps;++step){
        previous=x;const double dt=adaptiveFlowStep(field,x,h,o.flow,std::numeric_limits<double>::infinity());
        const auto u0=minus({previous[0],previous[1]},s.centre),u1=minus({x[0],x[1]},s.centre);
        // A step may not sweep more than 60 degrees about the centre.
        if(dot(u0,u1)<=.5*std::sqrt(dot(u0,u0)*dot(u1,u1))){x=previous;h=dt*.5;continue;}
        winding+=std::atan2(cross(u0,u1),dot(u0,u1));const double before=t;t+=dt;
        const double g0=g(previous),g1=g(x);
        if(!(g0<0&&g1>=0&&dot(d,u1)>0&&std::abs(winding)>kPi))continue;
        // Safeguarded Newton for the section crossing inside the accepted step.
        double lo=0,hi=dt,tau=dt*(-g0)/(g1-g0);V k(3);trial=previous;
        for(int it=0;it<80;++it){
            extrapolatedStep(field,previous,tau,trial);const double value=g(trial);
            if(value<0)lo=tau;else hi=tau;
            field(trial,k);const double slope=orient*(n[0]*k[0]+n[1]*k[1]);
            double next=slope!=0?tau-value/slope:.5*(lo+hi);if(!(next>lo&&next<hi))next=.5*(lo+hi);
            const bool done=std::abs(next-tau)<=4*std::numeric_limits<double>::epsilon()*(before+dt);tau=next;if(done)break;
        }
        extrapolatedStep(field,previous,tau,trial);
        const double closure=std::hypot(trial[0]-start[0],trial[1]-start[1]);
        if(closure>o.closureTolerance*radius)continue;
        PeriodicOrbit r;r.energy=energy;r.period=before+tau;r.frequency=2*kPi/r.period;
        r.action=std::abs(trial[2])/(2*kPi);r.start=start;r.orientation=trial[2]>=0?1:-1;
        r.closureError=closure;r.energyDrift=std::abs(s.hamiltonian({trial[0],trial[1]})-energy);r.steps=step+1;return r;
    }
    throw std::runtime_error("Orbit did not close within the step budget");
}
std::vector<PlanarPoint> orbitSamples(const PlanarOrbitSystem& s,const PeriodicOrbit& orbit,const OrbitFamilyOptions& o){
    validate(s);if(o.phaseSamples<3||o.phaseSamples>1000000)throw std::invalid_argument("Orbit needs 3..1000000 phase samples");
    const auto field=planar(s,1);std::vector<PlanarPoint> out;out.reserve(o.phaseSamples);V x{orbit.start[0],orbit.start[1]};
    const double step=orbit.period/double(o.phaseSamples);
    for(std::size_t j=0;j<o.phaseSamples;++j){out.push_back({x[0],x[1]});if(j+1<o.phaseSamples)integrateFlow(field,x,step,o.flow);}
    return out;
}
OrbitFourier orbitFourier(const std::vector<PlanarPoint>& samples,double frequency,int maximumHarmonic){
    std::vector<std::complex<double>> a,b;a.reserve(samples.size());b.reserve(samples.size());
    for(const auto& p:samples){a.emplace_back(p[0]);b.emplace_back(p[1]);}
    return {periodicFourierSeries(a,frequency,maximumHarmonic),periodicFourierSeries(b,frequency,maximumHarmonic)};
}
double orbitPhase(const PlanarOrbitSystem& s,const PlanarPoint& point,const OrbitFamilyOptions& o){
    const auto orbit=periodicOrbit(s,s.hamiltonian(point),o);auto grid=o;grid.phaseSamples=std::max<std::size_t>(o.phaseSamples,64);
    const auto samples=orbitSamples(s,orbit,grid);std::size_t best=0;double nearest=std::numeric_limits<double>::infinity();
    for(std::size_t j=0;j<samples.size();++j){const auto u=minus(samples[j],point);const double r=dot(u,u);if(r<nearest){nearest=r;best=j;}}
    const auto forward=planar(s,1),backward=planar(s,-1);
    const V base{samples[best][0],samples[best][1]};double offset=0;V y=base;
    // Gauss-Newton on the arc-time offset from the nearest sample.
    for(int it=0;it<24;++it){
        y=base;if(offset>0)integrateFlow(forward,y,offset,o.flow);else if(offset<0)integrateFlow(backward,y,-offset,o.flow);
        const auto f=s.drift({y[0],y[1]});const double speed=dot(f,f);if(!(speed>0))break;
        const double delta=-dot(f,minus({y[0],y[1]},point))/speed;offset+=delta;
        if(std::abs(delta)<1e-13*orbit.period)break;
    }
    const double phase=2*kPi*(double(best)/double(samples.size())+offset/orbit.period);
    const double wrapped=std::fmod(phase,2*kPi);return wrapped<0?wrapped+2*kPi:wrapped;
}
OrbitFamilyTable orbitFamily(const PlanarOrbitSystem& s,const V& energies,const OrbitFamilyOptions& o){
    requireIncreasing(energies,"Orbit energies");OrbitFamilyTable t;t.energies=energies;
    const double lo=energies.front(),hi=energies.back(),h=o.derivativeStep*(hi-lo);
    if(!(h>0))throw std::invalid_argument("Derivative step must be positive");
    const auto omega=[&](double e){return periodicOrbit(s,e,o).frequency;};
    for(double e:energies){
        auto orbit=periodicOrbit(s,e,o);const double w=orbit.frequency;double slope;
        if(e-2*h>=lo&&e+2*h<=hi)slope=(-omega(e+2*h)+8*omega(e+h)-8*omega(e-h)+omega(e-2*h))/(12*h);
        else if(e+2*h<=hi)slope=(-3*w+4*omega(e+h)-omega(e+2*h))/(2*h);
        else slope=(3*w-4*omega(e-h)+omega(e-2*h))/(2*h);
        t.periods.push_back(orbit.period);t.frequencies.push_back(w);t.actions.push_back(orbit.action);
        t.frequencySlopes.push_back(slope);t.orbits.push_back(orbit);
    }
    return t;
}
FrequencyBranchReport frequencyBranches(const V& e,const V& w){
    requireIncreasing(e,"Branch energies");if(e.size()!=w.size())throw std::invalid_argument("Branch table sizes differ");
    FrequencyBranchReport r;std::size_t start=0;int direction=0;
    const auto close=[&](std::size_t last){if(direction==0)return;
        r.branches.push_back({start,last,direction,std::min(w[start],w[last]),std::max(w[start],w[last])});};
    for(std::size_t i=0;i+1<e.size();++i){
        if(!std::isfinite(w[i])||!std::isfinite(w[i+1]))throw std::invalid_argument("Branch frequencies must be finite");
        const int sign=w[i+1]>w[i]?1:(w[i+1]<w[i]?-1:0);
        if(sign==0)++r.flatSegments;
        if(sign!=direction){close(i);start=i;direction=sign;}
    }
    close(e.size()-1);r.monotone=r.branches.size()==1&&r.flatSegments==0;return r;
}
std::vector<FrequencyRoot> invertFrequency(const V& e,const V& w,double target,const FrequencyEvaluator& evaluator){
    if(!std::isfinite(target))throw std::invalid_argument("Target frequency must be finite");
    const auto report=frequencyBranches(e,w);const auto slopes=tableDerivative(e,w);std::vector<FrequencyRoot> out;
    for(std::size_t b=0;b<report.branches.size();++b){
        const auto& branch=report.branches[b];if(target<branch.lower||target>branch.upper)continue;
        for(std::size_t i=branch.first;i<branch.last;++i){
            const double a0=w[i]-target,a1=w[i+1]-target;
            // A root on a shared turning node belongs to the earlier branch only.
            if(a0==0&&i==branch.first&&b>0&&report.branches[b-1].last==i)continue;
            const bool last=i+1==branch.last;
            if(!((a0<=0&&a1>0)||(a0>=0&&a1<0)||(last&&a1==0)))continue;
            FrequencyRoot r;r.branch=b;
            if(a0==0)r.energy=e[i];else if(a1==0)r.energy=e[i+1];
            else if(evaluator){
                double lo=e[i],hi=e[i+1],flo=a0,fhi=a1,x=lo;
                for(int it=0;it<60;++it){
                    x=(lo*fhi-hi*flo)/(fhi-flo);const double fx=evaluator(x)-target;
                    if(fx==0||hi-lo<=1e-13*(1+std::abs(x)))break;
                    if((fx<0)==(flo<0)){lo=x;flo=fx;fhi*=.5;}else{hi=x;fhi=fx;flo*=.5;}
                }
                r.energy=x;r.refined=true;
            }else r.energy=e[i]-a0*(e[i+1]-e[i])/(a1-a0);
            r.frequency=evaluator?evaluator(r.energy):target;
            const double fraction=(r.energy-e[i])/(e[i+1]-e[i]);
            r.slope=slopes[i]+fraction*(slopes[i+1]-slopes[i]);
            if(evaluator){const double h=std::min({r.energy-e.front(),e.back()-r.energy,1e-5*(e.back()-e.front())});
                if(h>0)r.slope=(evaluator(r.energy+h)-evaluator(r.energy-h))/(2*h);}
            out.push_back(r);break;
        }
    }
    return out;
}
std::vector<ShellRoot> shellLadder(const V& e,const V& w,double drive,int maximum,const FrequencyEvaluator& evaluator){
    if(!(drive>0)||!std::isfinite(drive))throw std::invalid_argument("Drive frequency must be positive");
    if(maximum<1||maximum>100000)throw std::invalid_argument("Shell harmonic limit out of range");
    std::vector<ShellRoot> out;
    for(int n=1;n<=maximum;++n)for(const auto& r:invertFrequency(e,w,drive/n,evaluator))out.push_back({n,r});
    return out;
}
SeparatrixFit fitSeparatrixPeriod(const V& e,const V& T,double separatrix){
    if(e.size()!=T.size()||e.size()<3)throw std::invalid_argument("Separatrix fit needs at least three matching samples");
    double sx=0,sy=0,sxx=0,sxy=0;const double n=double(e.size());
    for(std::size_t i=0;i<e.size();++i){
        const double delta=std::abs(separatrix-e[i]);
        if(!(delta>0)||!std::isfinite(T[i])||!(T[i]>0))throw std::invalid_argument("Separatrix fit needs positive periods off the separatrix");
        const double x=std::log(delta);sx+=x;sy+=T[i];sxx+=x*x;sxy+=x*T[i];
    }
    const double det=n*sxx-sx*sx;if(!(std::abs(det)>0))throw std::invalid_argument("Separatrix fit needs distinct energy gaps");
    const double slope=(n*sxy-sx*sy)/det,intercept=(sy-slope*sx)/n;
    if(!(slope<0))throw std::invalid_argument("Period does not grow toward the separatrix");
    SeparatrixFit fit;fit.lambda=-1/slope;fit.constant=std::exp(intercept*fit.lambda);fit.samples=e.size();
    double rss=0;for(std::size_t i=0;i<e.size();++i){const double r=T[i]-(intercept+slope*std::log(std::abs(separatrix-e[i])));rss+=r*r;}
    fit.rmsResidual=std::sqrt(rss/n);return fit;
}
}
