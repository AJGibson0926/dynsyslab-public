#include <dynsyslab/core/TorusAnalysis.hpp>
#include <cmath>
#include <stdexcept>

// Action-angle product Haar sampling and invariance equation DK(theta) omega=f(K).
// Arnold, Mathematical Methods of Classical Mechanics, chapter 10.
// Numerical residual only; rigorous bounds require a-posteriori KAM hypotheses.
namespace dynsyslab::core {
TorusSample sampleInvariantTorus(const ODESystem&system,const AnalyticTrajectory&trajectory,
    const std::vector<Observable>&obs,const std::vector<std::size_t>&shape,double t0){
    if(!trajectory||shape.empty()||shape.size()>8||!std::isfinite(t0))throw std::invalid_argument("Valid analytic torus and mesh required");
    if(system.configuration().driver.amplitude!=0)throw std::invalid_argument("Autonomous invariant-torus sampling requires zero prescribed driver");
    auto map=trajectory.torus;auto w=trajectory.angularFrequencies;const double tau=2*std::acos(-1.);
    if(!map&&std::isfinite(trajectory.period)&&trajectory.period>0){w={tau/trajectory.period};map=[=](const std::vector<double>&theta){return trajectory.evaluate(t0+theta[0]/w[0]);};}
    if(!map||w.size()!=shape.size())throw std::invalid_argument("No matching invariant-torus parametrization/period is available");
    for(double v:w)if(!std::isfinite(v))throw std::invalid_argument("Nonfinite torus frequency");
    const auto d=system.descriptor().stateDimension;std::size_t count=1;
    for(auto n:shape){if(n<3||n>1000000||count>1000000/n)throw std::invalid_argument("Torus product mesh exceeds limits");count*=n;}
    if(!d||d>10000000/count||obs.empty()||obs.size()>10000000/count)throw std::invalid_argument("Torus sample storage exceeds limits");
    for(const auto&o:obs)if(o.stateDimension()!=d||o.outputDimension()!=1)throw std::invalid_argument("Torus sampling expects scalar compatible observables");
    TorusSample out{shape,w,{}, {},d,obs.size()};out.states.reserve(count*d);out.observables.reserve(count*obs.size());
    auto state=[&](const std::vector<double>&theta){auto x=map(theta);if(x.size()!=d||system.invalid(t0,x.data()))throw std::invalid_argument("Torus map left the system domain");
        for(auto v:x)if(!std::isfinite(v))throw std::invalid_argument("Nonfinite torus state");return x;};
    for(std::size_t index=0;index<count;++index){std::vector<double> theta(shape.size());auto remaining=index;
        for(std::size_t a=shape.size();a-->0;){theta[a]=tau*double(remaining%shape[a])/double(shape[a]);remaining/=shape[a];}
        const auto x=state(theta);out.states.insert(out.states.end(),x.begin(),x.end());
        for(const auto&o:obs){auto v=o.evaluate(x.data())[0];if(!std::isfinite(v.real())||!std::isfinite(v.imag()))throw std::invalid_argument("Nonfinite torus observable");out.observables.push_back(v);}
        // Fourth-order directional difference in time, with bounded angular steps.
        double speed=0;for(auto v:w)speed=std::max(speed,std::abs(v));const double h=1e-4/std::max(1.,speed);
        auto moved=[&](double dt){auto phi=theta;for(std::size_t a=0;a<w.size();++a)phi[a]+=dt*w[a];return state(phi);};
        const auto mm=moved(-2*h),m=moved(-h),p=moved(h),pp=moved(2*h);std::vector<double> f(d);system.rhsWithInput(t0,x.data(),0,f.data());
        double error=0,norm=0;for(std::size_t j=0;j<d;++j){error=std::hypot(error,(-pp[j]+8*p[j]-8*m[j]+mm[j])/(12*h)-f[j]);norm=std::hypot(norm,f[j]);}
        out.maximumInvarianceResidual=std::max(out.maximumInvarianceResidual,error/(1+norm));
        for(std::size_t a=0;a<w.size();++a){auto phi=theta;phi[a]+=tau;const auto y=state(phi);double err=0,size=0;
            for(std::size_t j=0;j<d;++j){err=std::hypot(err,y[j]-x[j]);size=std::hypot(size,x[j]);}
            out.maximumWrapResidual=std::max(out.maximumWrapResidual,err/(1+size));}
    }return out;
}
}
