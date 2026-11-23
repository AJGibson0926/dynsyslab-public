#include <dynsyslab/core/PeriodicTransport.hpp>
#include <dynsyslab/core/Integrator.hpp>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace dynsyslab;
namespace {
void near(double a,double b,double tol=1e-8){if(!std::isfinite(a)||std::abs(a-b)>tol*(1+std::abs(b)))throw std::runtime_error("Periodic transport mismatch: "+std::to_string(a)+" vs "+std::to_string(b));}
template<class F>void rejects(F f){bool caught=false;try{f();}catch(const std::exception&){caught=true;}if(!caught)throw std::runtime_error("Expected rejection");}
math::PeriodicAtlas polar(){std::vector<math::PeriodicAtlasKnot> knots;
    for(int i=0;i<=64;++i){const double a=1.+i/64.;
        knots.push_back({a,{{1,{{-1,a/2},{1,a/2}}},{1,{{-1,{0,-a/2}},{1,{0,a/2}}}}}});
    }return math::PeriodicAtlas(knots);
}
}
int main(){try{
    const auto atlas=polar();const math::PeriodicChart chart(atlas,256);
    const double a=1.493,phase=.731,h=1e-6;const auto v=chart.evaluate(a,phase);
    near(v.value[0],a*std::cos(phase),2e-9);near(v.value[1],-a*std::sin(phase),2e-9);
    near(v.parameterDerivative[0],std::cos(phase),2e-9);
    for(int c=0;c<2;++c){
        near(v.parameterDerivative[c],(chart.evaluate(a+h,phase).value[c]-chart.evaluate(a-h,phase).value[c])/(2*h));
        near(v.phaseDerivative[c],(chart.evaluate(a,phase+h).value[c]-chart.evaluate(a,phase-h).value[c])/(2*h));
        near(chart.evaluate(a,2*std::acos(-1.)).value[c],chart.evaluate(a,0).value[c]);
    }
    const auto b=math::pullbackPlanarField(v,{0,1},{0,1});
    near(b.coordinates[0],-std::sin(phase),1e-6);near(b.coordinates[1],-std::cos(phase)/a,1e-6);
    near(b.normalizedDeterminant,-1,1e-9);
    rejects([&]{chart.evaluate(.9,0);});rejects([&]{math::pullbackPlanarField(v,{0,0},{0,1});});
    auto bad=v;bad.phaseDerivative=bad.parameterDerivative;rejects([&]{math::pullbackPlanarField(bad,{0,1},{0,1});});
    core::PeriodicTransportOptions o;o.phaseSamples=1024;o.inputMaximumHarmonic=127;
    int calls=0;core::PeriodicTransportSystem model(atlas,[&](auto){++calls;return std::array<double,2>{0,1};},o);
    const int buildCalls=calls;core::IntegrationOptions solver;solver.finalTime=8;solver.outputStep=.01;solver.maximumStep=.02;solver.relativeTolerance=1e-10;solver.absoluteTolerance=1e-12;
    auto folded=atlas.knots();for(auto& knot:folded){const double r=2+std::pow(knot.parameter-1.37,2);
        for(auto& channel:knot.channels)for(auto& [n,c]:channel.coefficients)c*=r/knot.parameter;}
    rejects([&]{core::PeriodicTransportSystem badMap(math::PeriodicAtlas(folded),[](auto){return std::array<double,2>{0,1};},o);});
    const auto x0=model.physicalState(std::array<double,2>{a,phase}.data());
    for(double damping:{0.,.02}){
        auto controlled=core::integrate(model,{a,phase},1,solver,[&](std::size_t,double,const double* z){return .03-damping*model.physicalState(z)[1];});
        const double w=std::sqrt(1-damping*damping/4),c=x0[0]-.03,d=(x0[1]+damping*c/2)/w;
        for(std::size_t i=0;i<controlled.sampleCount;++i){const double t=controlled.time[i],decay=std::exp(-damping*t/2);
            const double q=.03+decay*(c*std::cos(w*t)+d*std::sin(w*t));
            const double p=decay*((-damping*c/2+w*d)*std::cos(w*t)+(-damping*d/2-w*c)*std::sin(w*t));
            const auto x=model.physicalState(controlled.state(0,i));near(x[0],q,2e-6);near(x[1],p,2e-6);
        }
    }
    if(calls!=buildCalls)throw std::runtime_error("Physical field used during reduced rollout");
    std::cout<<"Periodic chart derivatives, input transport, forced and feedback solutions passed\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
