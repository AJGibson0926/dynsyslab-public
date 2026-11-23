#include <dynsyslab/core/System.hpp>
#include <dynsyslab/core/CoefficientShift.hpp>
#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>

using namespace dynsyslab::core;
void check(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
void verify(const std::string& name,std::vector<double> initial={}){
    SystemConfiguration c;c.identifier=name;c.initialState=initial;auto system=createSystem(c);
    const auto x0=system->initialState();const auto families=system->analyticSolutions();
    check(!families.empty(),"Missing analytic family");auto a=families[0].instantiate(x0,1.3);
    if(!a)throw std::runtime_error(name+": "+a.reason);
    const auto H=system->observableCatalog().at("total_energy");
    const double h0=H.evaluate(x0.data())[0].real();
    for(int k=0;k<60;++k){const double t=1.3+.137*k,step=1e-5;
        auto x=a.evaluate(t),xm=a.evaluate(t-step),xp=a.evaluate(t+step);std::vector<double> f(x.size());
        system->rhsWithInput(t,x.data(),0,f.data());
        for(std::size_t j=0;j<x.size();++j){
            if(std::abs((xp[j]-xm[j])/(2*step)-f[j])>1e-6*(1+std::abs(f[j])))throw std::runtime_error(name+": ODE residual");
            if(k==0)check(std::abs(x[j]-x0[j])<2e-10,"IC mismatch");}
        check(std::abs(H.evaluate(x.data())[0].real()-h0)<2e-9*(1+std::abs(h0)),"Energy changed");
    }
    if(a.period>0){const auto x=a.evaluate(1.3+a.period);for(std::size_t j=0;j<x.size();++j)check(std::abs(x[j]-x0[j])<1e-9,"Period mismatch");}
    if(a.period>0) {
        const auto covariance=checkPeriodicCoefficientShift(families[0],x0,
            {system->observableCatalog().at("state")},512,31,{.13*a.period,-.37*a.period},1.3);
        for(const auto& row:covariance.shifts) {
            check(row.relativeFrequencyError<1e-10&&row.initialStateError<1e-9,"Shifted IC or invariant period changed");
            for(std::size_t j=0;j<row.coefficientL2Error.size();++j)
                check(row.coefficientL2Error[j]<1e-7*(1+row.referenceL2Norm[j]),"Fourier coefficient time-shift identity");
        }
    }
    if(a.torus){std::vector<double> angles=a.angularFrequencies;for(auto&v:angles)v*=.17;
        const auto x=a.torus(angles),y=a.evaluate(1.47);for(std::size_t j=0;j<x.size();++j)check(std::abs(x[j]-y[j])<1e-9,"Torus flow mismatch");}
    c.parameters["damping"]=.01;auto damped=createSystem(c);check(!damped->analyticSolutions()[0].instantiate(x0,0),"Damped orbit accepted as invariant torus");
    std::cout<<name<<" passed\n";
}
int main(){try{
    for(auto name:{"morse_oscillator","isotonic_oscillator","poschl_teller","cosh_oscillator","radial_twist","kepler_planar","euler_top","coupled_duffing_hamiltonian"})verify(name);
    verify("euler_top",{1,.3,.2});verify("euler_top",{-.2,.3,-1});verify("euler_top",{-1,.3,.2});
    verify("morse_oscillator",{-.3,.2});verify("cosh_oscillator",{.3,.4});verify("kepler_planar",{1,.3,-.2,-.7});
    verify("kepler_planar",{1,0,0,1});
    {
        SystemConfiguration config;config.identifier="radial_twist";config.parameters["omega"]=1;config.parameters["twist"]=1;
        const auto system=createSystem(config);const auto catalog=system->observableCatalog();
        for(double I:{.5,1.,3.,6.5}){
            const std::vector<double> initial{std::sqrt(2*I),0};
            auto solution=system->analyticSolutions()[0].instantiate(initial,0);
            check(std::abs(catalog.at("action").evaluateReal(initial.data())[0]-I)<1e-13,"Radial action");
            check(std::abs(catalog.at("hamiltonian_frequency").evaluateReal(initial.data())[0]-(1+I))<1e-13,"Radial frequency");
            auto phase=[&](double t){const auto x=solution.evaluate(t);return std::complex<double>{
                catalog.at("phase_cosine").evaluateReal(x.data())[0],catalog.at("phase_sine").evaluateReal(x.data())[0]};};
            const auto before=phase(-1),now=phase(0),after=phase(1);
            check(std::abs(after*before-now*now)<1e-12,"Exact two-sample phase recurrence");
            const auto x=solution.evaluate(.37);const auto z=phase(.37);
            const std::complex<double> dz{catalog.at("phase_cosine_dot").evaluateReal(x.data())[0],
                catalog.at("phase_sine_dot").evaluateReal(x.data())[0]};
            check(std::abs(dz-std::complex<double>(0,1+I)*z)<1e-12,"Phase derivative links");
        }
    }
    SystemConfiguration c;c.identifier="coupled_duffing_hamiltonian";c.parameters["coupling"]=.1;
    auto s=createSystem(c);check(!s->analyticSolutions()[0].instantiate(s->initialState(),0),"Perturbed system accepted as separable");
    AnalyticSolutionFamily broken{{},[](const AnalyticState&,double t0) {
        return AnalyticTrajectory{[=](double t){return AnalyticState{std::cos(t-t0),std::sin(t-t0)};},{},2*std::acos(-1.)};
    }};
    const auto bad=checkPeriodicCoefficientShift(broken,{1,0},{Observable::identity(2,{"q","p"},{"q","p"})},64,3,{.7});
    check(bad.shifts[0].coefficientL2Error[0]>.1&&bad.shifts[0].initialStateError>.1,
        "Coefficient check detects a family that discards initial phase");
    std::cout<<"Integrable system tests passed\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
