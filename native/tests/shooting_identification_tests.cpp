#include <dynsyslab/core/ShootingIdentification.hpp>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace dynsyslab;
namespace {
void require(bool b,const char* message){if(!b)throw std::runtime_error(message);}
template<class F>void rejects(F f){bool caught=false;try{f();}catch(const std::exception&){caught=true;}require(caught,"Expected rejection");}
}
int main(){try{
    core::SystemConfiguration config;config.identifier="harmonic_oscillator";auto system=core::createSystem(config);
    const double x0=1.31,v0=-.27,b=.2,feedthrough=.7;
    std::vector<double> time,y;
    for(int j=0;j<=100;++j){const double t=j*.03;time.push_back(4+t);y.push_back(b+(x0-b)*std::cos(t)+v0*std::sin(t)+feedthrough*b);}
    core::IntegrationOptions solver;solver.relativeTolerance=1e-11;solver.absoluteTolerance=1e-13;solver.maximumStep=.03;
    const auto initial=[](const std::vector<double>& p){return p;};
    const auto read=[&](double,const double* x,double u){return std::vector<double>{x[0]+feedthrough*u};};
    const auto input=[&](double){return b;};
    std::vector<core::ShootingStart> starts{{{1.,.1},{.5,-1},{2,1}},{{1.7,-.5},{.5,-1},{2,1}}};
    auto fit=core::identifyByShooting(*system,time,y,{.01},starts,initial,read,input,solver);
    require(std::abs(fit.best.parameters[0]-x0)<1e-8&&std::abs(fit.best.parameters[1]-v0)<1e-8,"Exact driven oscillator IC mismatch");
    require(fit.best.squaredResidual<1e-14,"Exact driven oscillator residual mismatch");
    // A stale failed start cannot suppress a good candidate or appear successful.
    auto selective=[&](const std::vector<double>& p){if(p[0]>1.6)throw std::domain_error("unavailable initial chart");return p;};
    fit=core::identifyByShooting(*system,time,y,{.01},starts,selective,read,input,solver);
    require(fit.attempts[0].valid&&!fit.attempts[1].valid&&!fit.attempts[1].failure.empty(),"Failed candidate diagnostics lost");
    const auto wrong=[](double,const double* x,double){return std::vector<double>{x[0]};};
    fit=core::identifyByShooting(*system,time,y,{.01},starts,initial,wrong,input,solver);
    require(fit.best.squaredResidual>.1,"Direct input observation term was not tested");
    auto bad=time;bad[20]+=.01;
    rejects([&]{core::identifyByShooting(*system,bad,y,{.01},starts,initial,read,input,solver);});
    rejects([&]{core::identifyByShooting(*system,time,y,{0},starts,initial,read,input,solver);});
    std::cout<<"Controlled shooting: exact solution, input observation, failed starts and validation passed\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
