#include <dynsyslab/core/System.hpp>
#include <dynsyslab/core/Integrator.hpp>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace dynsyslab::core;
namespace {
void check(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
void near(double a,double b,double tolerance=1e-10){
    check(std::isfinite(a)&&std::isfinite(b)&&std::abs(a-b)<tolerance*(1+std::abs(b)),"Multiwell numeric mismatch");
}
void verify(const std::string& id,const std::string& preset="") {
    SystemConfiguration config;config.identifier=id;config.preset=preset;
    auto system=createSystem(config);
    auto catalog=system->observableCatalog();
    const auto& energy=catalog.at("total_energy");
    for(double a:{-.8,-.2,0.,.4,1.1}) for(double b:{-.3,0.,.2}) {
        const double x[]{a,b};double f[2],controlled[2];
        system->rhsWithInput(0,x,0,f);system->rhsWithInput(0,x,.37,controlled);
        const auto gradient=energy.jacobian(x);
        near((gradient[0]*f[0]+gradient[1]*f[1]).real(),catalog.at("energy_dot").evaluateReal(x)[0]);
        near((gradient[0]*(controlled[0]-f[0])+gradient[1]*(controlled[1]-f[1])).real(),
             .37*catalog.at("input_power_gain").evaluateReal(x)[0]);
        auto jac=system->jacobian(0,x);
        for(size_t j=0;j<2;++j){
            double xp[]{a,b},xm[]{a,b},fp[2],fm[2];xp[j]+=1e-5;xm[j]-=1e-5;
            system->rhsWithInput(0,xp,0,fp);system->rhsWithInput(0,xm,0,fm);
            for(size_t i=0;i<2;++i)near(jac[2*i+j],(fp[i]-fm[i])/2e-5,2e-9);
        }
        for(const auto& link:catalog.timeDerivatives()){
            const auto actual=catalog.at(link.first).flowDerivative(system->driftExpressions()).evaluateReal(x);
            const auto expected=catalog.at(link.second.target).evaluateReal(x);
            near(actual.at(0),expected.at(0));
        }
    }
    IntegrationOptions options;options.finalTime=5;options.outputStep=.05;
    options.relativeTolerance=1e-11;options.absoluteTolerance=1e-13;
    const auto initial=system->initialState();
    const double initialEnergy=energy.evaluateReal(initial.data())[0];
    auto solution=integrate(*system,initial,1,options);
    for(size_t n=0;n<solution.sampleCount;++n)
        near(energy.evaluateReal(solution.state(0,n))[0],initialEnergy,2e-9);
}
}
int main(){try{
    for(const auto& p:{"conservative_chain","left_libration","right_libration","exterior_orbit","equal_depth_wells"})
        verify("triple_well_duffing",p);
    for(const auto& p:{"right_leaf","upper_left_leaf","lower_left_leaf","exterior_orbit","separatrix"})
        verify("three_leaf_hamiltonian_clover",p);
    SystemConfiguration c;c.identifier="three_leaf_hamiltonian_clover";
    auto clover=createSystem(c);double f[2];const double origin[]{0,0},center[]{.75,0},tip[]{1,0};
    clover->rhsWithInput(0,origin,0,f);near(f[0],0);near(f[1],0);
    for(double v:clover->jacobian(0,origin))near(v,0);
    clover->rhsWithInput(0,center,0,f);near(f[0],0);near(f[1],0);
    clover->rhsWithInput(0,tip,0,f);near(f[0],0);near(f[1],-1);
    near(clover->observableCatalog().at("total_energy").evaluateReal(center)[0],-27./256);
    c.identifier="triple_well_duffing";c.parameters["damping"]=.2;
    auto damped=createSystem(c);const double x[]{.3,.4};
    near(damped->observableCatalog().at("energy_dot").evaluateReal(x)[0],-.2*.4*.4);
    c.parameters["inner_saddle"]=1;
    bool rejected=false;try{(void)createSystem(c);}catch(const std::invalid_argument&){rejected=true;}
    check(rejected,"Invalid well geometry accepted");
    std::cout<<"Multiwell systems, observables, derivatives, energy and Jacobians passed\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
