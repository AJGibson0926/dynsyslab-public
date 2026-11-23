#include <dynsyslab/core/OnlineShooting.hpp>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace dynsyslab;
namespace {
void check(bool b,const char* message){if(!b)throw std::runtime_error(message);}
template<class F>void rejects(F f){try{f();}catch(const std::exception&){return;}throw std::runtime_error("Expected observer rejection");}
}
int main(){try{
    core::SystemConfiguration config;config.identifier="harmonic_oscillator";auto model=core::createSystem(config),other=core::createSystem(config);
    core::OnlineShootingOptions options;options.sampleTime=.03;options.firstFitIndex=20;options.fitEvery=10;options.maximumSamples=101;
    core::IntegrationOptions solver;solver.relativeTolerance=1e-11;solver.absoluteTolerance=1e-13;solver.maximumStep=.03;
    const std::vector<core::ShootingStart> starts{{{1,.1},{.5,-1},{2,1}}};
    const auto initialize=[](const auto& p){return p;};
    const auto read=[](double,const double* x,double u){return std::vector<double>{x[0]+.7*u};};
    int providerCalls=0;
    const auto provider=[&](const auto& t,const auto& y){++providerCalls;check(t.size()==21&&y.size()==21,"Initializer saw a future prefix");return starts;};
    core::OnlineShootingObserver observer(*model,{.01},{},initialize,read,options,solver,{},provider);
    core::OnlineShootingObserver perturbed(*other,{.01},starts,initialize,read,options,solver);
    rejects([&]{observer.forecast(1,0,10);});
    double x=1.31,v=-.27,u=0,maxError=0,lateDifference=0;
    for(std::size_t k=0;k<101;++k){
        const double time=k*.03;
        const auto sample=observer.observe(time,{x+.7*u},u);
        const auto changed=perturbed.observe(time,{x+.7*u+(k>=60?.02:0)},u);
        check(sample.ready==(k>=20),"Acquisition readiness mismatch");
        if(k<20)check(sample.state.empty()&&!sample.informationTime,"Pre-acquisition truth leaked");
        else{
            check(*sample.informationTime<=time,"Observer used future measurements");
            maxError=std::max(maxError,std::hypot(sample.state[0]-x,sample.state[1]-v));
            if(k<60)check(sample.state==changed.state,"Later sensor data changed an earlier estimate");
            else lateDifference=std::max(lateDifference,std::abs(sample.state[0]-changed.state[0]));
            if(sample.fit)check(sample.fit->best.squaredResidual<1e-12,"Exact switched sensor fit mismatch");
        }
        const double next=k<20?.1:std::clamp(-.2*sample.state[1],-.5,.5);
        if(sample.ready){const auto forecast=observer.forecast(time+.03,next,8);const auto last=forecast.state(0,forecast.sampleCount-1);
            check(std::abs(last[0]-(next+(x-next)*std::cos(.03)+v*std::sin(.03)))<1e-7,"Causal forecast is not the exact held-input solution");}
        const double nx=next+(x-next)*std::cos(.03)+v*std::sin(.03);v=-(x-next)*std::sin(.03)+v*std::cos(.03);x=nx;u=next;
    }
    check(providerCalls==1&&maxError<1e-7&&lateDifference>1e-4,"Online acquisition/recovery/future sensitivity mismatch");
    rejects([&]{observer.observe(3.03,{x},u);});
    core::OnlineShootingObserver clock(*model,{.01},starts,initialize,read,options,solver);
    (void)clock.observe(2,{1},0);rejects([&]{clock.observe(2.031,{1},0);});
    std::cout<<"Online shooting: exact closed-loop plant, sensor timing, acquisition, causal forecasts and future independence passed\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
