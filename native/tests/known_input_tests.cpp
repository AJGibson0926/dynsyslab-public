#include <dynsyslab/core/KnownInput.hpp>
#include <dynsyslab/core/InputAffineSensor.hpp>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace dynsyslab::core;
namespace {
void check(bool b,const char* m){if(!b)throw std::runtime_error(m);}
template<class F>void rejects(F f){try{f();}catch(const std::exception&){return;}throw std::runtime_error("Expected input rejection");}
}
int main(){try{
    const std::vector<double> time{2,2.173,2.1731,2.8,3.0},values{.1,200,-.3,.4,7};
    auto input=KnownScalarInput::recorded(time,values,-.2);
    check(input(2)==.1&&input.measurement(2)==-.2,"Initial input association");
    check(input(2.173)==200&&input.measurement(2.173)==.1,"Switch input association");
    check(input(3)==7&&input.measurement(3)==.4,"Final sample association");
    auto after=KnownScalarInput::recorded(time,values,-.2,false);
    check(after.measurement(2.173)==200,"After-command association");
    rejects([&]{input(3.01);});rejects([&]{input(1.99);});
    rejects([&]{KnownScalarInput::recorded({0,0},{1,2},0);});
    SystemConfiguration config;config.identifier="harmonic_oscillator";auto system=createSystem(config);
    IntegrationOptions solver;solver.initialTime=2;solver.finalTime=3;solver.outputStep=.03;
    solver.relativeTolerance=1e-11;solver.absoluteTolerance=1e-13;
    const auto out=integrateKnownInput(*system,{1.2,-.4},solver,input);
    // Independent exact constant-forcing harmonic transitions, including a
    // narrow pulse far smaller than output_step. No numerical reference solve.
    for(std::size_t j=0;j<out.sampleCount;++j){double x=1.2,v=-.4,current=2;
        for(std::size_t k=0;k+1<time.size()&&current<out.time[j];++k){
            const double stop=std::min(out.time[j],time[k+1]),dt=stop-current,b=values[k];
            if(dt>0){const double nx=b+(x-b)*std::cos(dt)+v*std::sin(dt);v=-(x-b)*std::sin(dt)+v*std::cos(dt);x=nx;current=stop;}
        }
        check(std::abs(out.state(0,j)[0]-x)<2e-9&&std::abs(out.state(0,j)[1]-v)<2e-9,"Exact switched harmonic replay mismatch");
    }
    check(out.time.back()==3&&out.sampleCount==35,"Replay grid mismatch");
    solver.outputStep=.01;
    const auto roundoff=integrateKnownInput(*system,{1.2,-.4},solver,input);
    check(roundoff.sampleCount==101,"Decimal clock acquired a duplicate endpoint");
    for(std::size_t j=1;j<roundoff.time.size();++j)check(roundoff.time[j]>roundoff.time[j-1],"Replay clock not strictly increasing");
    solver.initialTime=2.2;solver.finalTime=2.7;
    const auto partial=integrateKnownInput(*system,{1,0},solver,input);
    check(std::abs(partial.state(0,partial.sampleCount-1)[0]-(-.3+1.3*std::cos(.5)))<1e-9,"Replay window inside a held interval");
    const auto catalog=system->observableCatalog();const double state[]{2,3};
    InputAffineSensor plain({catalog.at("position"),catalog.at("velocity")});
    InputAffineSensor scaled({catalog.at("position"),catalog.at("velocity")},{.5,-.2},catalog.at("velocity"));
    check(plain.evaluate(state,.4)==std::vector<double>({2,3}),"Default sensor changed");
    const auto y=scaled.evaluate(state,.4);
    check(std::abs(y[0]-2.6)<1e-14&&std::abs(y[1]-2.76)<1e-14,"Input-scaled sensor mismatch");
    rejects([&]{InputAffineSensor bad({catalog.at("position")},{1,2});});
    rejects([&]{InputAffineSensor bad({catalog.at("position")},{1},catalog.at("position"),1);});
    std::cout<<"Recorded ZOH replay, exact narrow-pulse integration, input timing and sensor tests passed\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
