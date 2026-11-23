#include <dynsyslab/core/SequentialShooting.hpp>
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
    const double x0=1.31,v0=-.27,b=.2;
    std::vector<double> time,y;
    for(int j=0;j<=300;++j){const double t=j*.03;time.push_back(4+t);y.push_back(b+(x0-b)*std::cos(t)+v0*std::sin(t)+.7*b);}
    const auto initialize=[](const std::vector<double>& p){return p;};
    const auto read=[](double,const double* x,double u){return std::vector<double>{x[0]+.7*u};};
    const auto input=[&](double){return b;};
    core::IntegrationOptions solver;solver.relativeTolerance=1e-11;solver.absoluteTolerance=1e-13;solver.maximumStep=.03;
    const std::vector<core::ShootingStart> starts{{{1.,.1},{.5,-1},{2,1}}};
    core::SequentialShootingOptions o;o.updateIndices={100,150,200};o.predictionSubsteps=4;
    const auto fit=core::sequentialShooting(*system,time,y,{.01},starts,initialize,read,input,o,solver);
    require(fit.updates.size()==3&&fit.time.size()==801,"Sequential dimensions mismatch");
    for(std::size_t j=0;j<fit.time.size();++j){const double t=fit.time[j]-4;
        require(fit.informationTime[j]<=fit.time[j],"Future information leaked");
        require(j==0||fit.time[j]>fit.time[j-1],"Duplicate output timestamp");
        require(std::abs(fit.states[j*2]-(b+(x0-b)*std::cos(t)+v0*std::sin(t)))<1e-7,"Exact sequential trajectory mismatch");
    }
    auto changed=y;for(std::size_t j=150;j<changed.size();++j)changed[j]+=.1;
    const auto perturbed=core::sequentialShooting(*system,time,changed,{.01},starts,initialize,read,input,o,solver);
    double laterDifference=0;
    for(std::size_t j=0;j<fit.time.size();++j){
        if(fit.time[j]<time[150])require(fit.states[j*2]==perturbed.states[j*2],"Future observation changed earlier prediction");
        else laterDifference=std::max(laterDifference,std::abs(fit.states[j*2]-perturbed.states[j*2]));
    }
    require(laterDifference>.001,"Later observations did not update predictions");
    auto bad=o;bad.updateIndices={0};rejects([&]{core::sequentialShooting(*system,time,y,{.01},starts,initialize,read,input,bad,solver);});
    bad.updateIndices={100,90};rejects([&]{core::sequentialShooting(*system,time,y,{.01},starts,initialize,read,input,bad,solver);});
    o.updateIndices={300};const auto terminal=core::sequentialShooting(*system,time,y,{.01},starts,initialize,read,input,o,solver);
    require(terminal.time.size()==1&&terminal.time.front()==time.back(),"Final-time update mismatch");
    std::cout<<"Sequential shooting: independent solution, future-data isolation, updates and clock boundaries passed\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
