// Analytic references for the RPE -> Hankel -> H-DMD -> KLQR chain. They replace the retired MATLAB
// exporters export_dynsyslab_rpe_hdmd_sindy.m and run_matlab_rpe_klqr_control_experiment.m, which
// compared the native chain with a legacy MATLAB suite; the fixed configuration below is the one recorded
// at the 2026-08-26 parity checkpoint (docs/notes/RPE_MATLAB_PARITY_CHECKPOINT_2026-08-26.md).
// Only stable public contracts are used: the system catalog, the integrator, buildHankelMatrix,
// the koopman-spectral query (model JSON passed back opaquely) and the numerical LQR solver.
//
// Rayleigh-Plesset: Lord Rayleigh, Phil. Mag. 34 (1917) 94-98; M. S. Plesset, J. Appl. Mech. 16 (1949)
//   277-282; linear oscillation and viscous damping: M. S. Plesset and A. Prosperetti, Annu. Rev. Fluid
//   Mech. 9 (1977) 145-185; C. E. Brennen, Cavitation and Bubble Dynamics (Oxford, 1995), Ch. 4.
// Hankel DMD: H. Arbabi and I. Mezic, SIAM J. Appl. Dyn. Syst. 16 (2017) 2096-2126; exact DMD:
//   J. H. Tu et al., J. Comput. Dyn. 1 (2014) 391-421; residuals: M. J. Colbrook and A. Townsend,
//   Comm. Pure Appl. Math. 77 (2024) 221-283.
// LQR: R. E. Kalman, Bol. Soc. Mat. Mexicana 5 (1960) 102-119; B. D. O. Anderson and J. B. Moore,
//   Optimal Control: Linear Quadratic Methods (Prentice Hall, 1990).
#include <dynsyslab/core/Hankel.hpp>
#include <dynsyslab/core/Integrator.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/System.hpp>
#include <dynsyslab/numerics/Control.hpp>
#include <dynsyslab/workflow/KoopmanSpectralQuery.hpp>

#include <cmath>
#include <complex>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
namespace core=dynsyslab::core;
namespace json=core::json;
using J=json::Value;
using cd=std::complex<double>;
constexpr double pi=3.14159265358979323846;

void require(bool condition,const std::string& message){if(!condition)throw std::runtime_error(message);}
void close(double actual,double expected,double tolerance,const std::string& label){
    if(!std::isfinite(actual)||std::abs(actual-expected)>tolerance)
        throw std::runtime_error(label+": actual="+std::to_string(actual)+" expected="+std::to_string(expected));
}
J array(const std::vector<double>& values){J::Array out;for(double v:values)out.emplace_back(v);return out;}
std::vector<double> numbers(const J& value){std::vector<double> out;for(const auto& v:value.array())out.push_back(v.number());return out;}

// The parity-checkpoint medium: air in water at 20 C (checkpoint table, SI units).
constexpr double R0=1e-5,kappa=4./3.,p0=1e5,pV=2330,rho=998,sigma=.0725,mu=.001,pG0=112170;

std::unique_ptr<core::ODESystem> rpe(double viscosity=mu){
    core::SystemConfiguration c;c.identifier="rpe";c.preset="air_water_20c";
    if(viscosity!=mu)c.parameters["mu_pa_s"]=viscosity;
    return core::createSystem(c);
}

// Natural-period scaling: tau = t f0 with the Minnaert frequency including surface tension,
// f0 = sqrt(3 kappa pG0 - 2 sigma/R0)/(2 pi R0 sqrt(rho)) (Brennen 1995, Sec. 4.2).
void checkpointConfiguration(){
    const auto system=rpe();const auto& c=system->configuration();
    const std::pair<const char*,double> parameters[]{{"equilibrium_radius_m",R0},{"ambient_pressure_pa",p0},
        {"rho_kg_m3",rho},{"mu_pa_s",mu},{"sigma_n_m",sigma},{"vapor_pressure_pa",pV},{"kappa",kappa}};
    for(const auto& [name,value]:parameters)
        if(c.parameters.count(name))close(c.parameters.at(name),value,1e-15*std::abs(value),std::string("checkpoint ")+name);
    const auto& g=c.dimensionlessGroups;
    const double stiffness=3*kappa*pG0-2*sigma/R0,f0=std::sqrt(stiffness)/(2*pi*R0*std::sqrt(rho));
    close(g.at("natural_frequency_hz"),f0,1e-9*f0,"Minnaert natural frequency");
    close(g.at("characteristic_time_s"),1/f0,1e-9/f0,"natural-period time scale");
    close(g.at("gas_pressure_ratio")*g.at("pressure_scale_pa"),pG0,1e-8*pG0,"equilibrium gas pressure pG0=p0-pV+2 sigma/R0");
    close(3*kappa*g.at("gas_pressure_ratio")-g.at("capillarity"),4*pi*pi,1e-10,"unit-period linear stiffness");
    close(g.at("viscosity"),4*mu/(rho*R0*R0*f0),1e-12,"viscous group 4 mu T/(rho R0^2)");
    require(system->descriptor().preferredCoordinates=="dimensionless","RPE solves in dimensionless coordinates");
}

core::NativeSolution solve(core::ODESystem& system,std::vector<double> x0,double final,double step,
                           core::IntegrationMethod method,double rtol,double atol){
    core::IntegrationOptions o;o.finalTime=final;o.outputStep=step;o.method=method;
    o.relativeTolerance=rtol;o.absoluteTolerance=atol;
    return core::integrate(system,x0,1,o);
}

// Canonical training trajectory of the checkpoint: r(0)=2.3333865986017348, BDF, rtol=atol=1e-9,
// tau in [0,25] at dtau=1e-3. The RPE first integral E (integrating factor r^3) obeys dE/dtau = -nu r v^2;
// with nu=0 it is conserved (Plesset 1949), otherwise its loss equals the viscous work.
void canonicalTrajectory(){
    const std::vector<double> x0{2.3333865986017348,0};
    for(double viscosity:{0.,mu}){
        auto system=rpe(viscosity);const auto catalog=system->observableCatalog();
        const auto s=solve(*system,x0,25,1e-3,core::IntegrationMethod::bdf,1e-9,1e-9);
        require(s.sampleCount==25001&&s.coordinateSystem=="dimensionless","checkpoint output grid");
        const auto energy=[&](std::size_t k){return catalog.at("mechanical_energy").evaluateReal(s.state(0,k))[0];};
        const double e0=energy(0);require(e0>0,"positive initial energy");
        double work=0,worst=0;
        for(std::size_t k=1;k<s.sampleCount;++k){
            const double a=catalog.at("energy_dissipation").evaluateReal(s.state(0,k-1))[0];
            const double b=catalog.at("energy_dissipation").evaluateReal(s.state(0,k))[0];
            work+=.5*(a+b)*(s.time[k]-s.time[k-1]);
            worst=std::max(worst,std::abs(energy(k)-e0-work)/e0);
        }
        require(worst<(viscosity==0?1e-6:2e-4),"RPE energy balance along the checkpoint trajectory ("+std::to_string(worst)+")");
        if(viscosity>0)require(energy(s.sampleCount-1)<.5*e0,"viscous ringdown loses energy");
    }
}

// Small-amplitude ringdown: r=1+delta with delta''+nu delta'+4 pi^2 delta=0 (Plesset & Prosperetti 1977),
// so delta=d0 e^{-nu tau/2}(cos w tau+(nu/2w) sin w tau), w=sqrt(4 pi^2-nu^2/4).
struct Linear { double rate,frequency; };
Linear linearRates(){const auto system=rpe();const double nu=system->configuration().dimensionlessGroups.at("viscosity");
    return {-nu/2,std::sqrt(4*pi*pi-nu*nu/4)};}

void linearRingdownAndHdmd(){
    constexpr double d0=1e-6,dtau=.01;
    const auto [rate,w]=linearRates();
    auto system=rpe();
    const auto s=solve(*system,{1+d0,0},6,dtau,core::IntegrationMethod::dormandPrince853,1e-13,1e-16);
    std::vector<double> radius;
    for(std::size_t k=0;k<s.sampleCount;++k){const double t=s.time[k];
        const double exact=d0*std::exp(rate*t)*(std::cos(w*t)-rate/w*std::sin(w*t));
        close(s.state(0,k)[0]-1,exact,5e-11,"linearized RPE ringdown");radius.push_back(s.state(0,k)[0]);}

    // Hankel matrix of the radius, as exported by the MATLAB wrapper: X[d,j]=r[j+d], Y one tick later.
    core::HankelWindow window;window.delays=20;window.columns=500;
    const auto h=core::buildHankelMatrix(radius,radius.size(),1,{{0,0,1.}},window,true);
    require(h.rows==20&&h.columns==500,"Hankel shape");
    for(std::size_t d:{std::size_t(0),std::size_t(7),std::size_t(19)})for(std::size_t j:{std::size_t(0),std::size_t(311),std::size_t(499)}){
        close(h.x[d*h.columns+j],radius[h.starts[j]+d],0,"Hankel X entry");
        close(h.y[d*h.columns+j],radius[h.starts[j]+d+1],0,"Hankel Y entry");}

    // H-DMD = DMD on Hankel snapshot pairs. r = 1 + damped cosine has exact rank three, with
    // discrete eigenvalues 1 and exp((rate +- i w) dtau).
    std::vector<double> x(h.rows*h.columns),y(x.size());
    for(std::size_t j=0;j<h.columns;++j)for(std::size_t d=0;d<h.rows;++d){
        x[j*h.rows+d]=h.x[d*h.columns+j];y[j*h.rows+d]=h.y[d*h.columns+j];}
    const cd mu1=std::exp(cd(rate,w)*dtau);
    J::Array grid{J::Array{1.,0.},J::Array{mu1.real(),mu1.imag()},J::Array{mu1.real(),-mu1.imag()},J::Array{.5,.5}};
    const auto result=dynsyslab::workflow::evaluateKoopmanSpectralQuery(J::Object{
        {"schema","dynsyslab.koopman-spectral-query/v1"},{"action","fit"},
        {"data",J::Object{{"dimension",double(h.rows)},{"x",array(x)},{"y",array(y)}}},
        {"settings",J::Object{{"method","edmd"},{"rank",3.},{"sample_interval",dtau},
            {"dictionary",J::Object{{"kind","linear"},{"constant",false}}},{"residual_grid",grid}}}});
    require(result.at("rank").number()==3,"H-DMD retained rank");
    const auto minimum=numbers(result.at("spectral").at("residuals").at("minimum"));
    require(minimum.size()==4,"residual grid");
    for(std::size_t i=0;i<3;++i)require(minimum[i]<1e-6,"analytic H-DMD eigenvalue "+std::to_string(i)+" residual "+std::to_string(minimum[i]));
    require(minimum[3]>1e-2,"a point off the spectrum keeps a large residual");
    // Free H-DMD forecast of the delay window from its first column only.
    const std::size_t steps=150;
    const auto rollout=numbers(dynsyslab::workflow::evaluateKoopmanSpectralQuery(J::Object{
        {"schema","dynsyslab.koopman-spectral-query/v1"},{"action","rollout"},{"model",result.at("model")},
        {"initial",array(std::vector<double>(x.begin(),x.begin()+h.rows))},{"steps",double(steps)}}).at("values"));
    for(std::size_t k=0;k<=steps;++k)close(rollout[k*h.rows],radius[h.starts[0]+k],1e-10,"H-DMD forecast of the radius");
}

// KLQR solves a continuous LQR on the lifted linear model each time its input map changes. The
// Riccati solver must give the closed form for the double integrator (Anderson & Moore 1990):
// Q=I, R=1: P=[[sqrt3,1],[1,sqrt3]], K=[1,sqrt3]; and a stabilizing gain on the linearized RPE with
// the pressure input map at the controlled equilibrium IC (1,0) used by the MATLAB KLQR experiment.
void riccati(){
    const auto lqr=dynsyslab::numerics::solveContinuousLqr({0,1,0,0},{0,1},{1,0,0,1},{1},2,1);
    const double s3=std::sqrt(3.);
    close(lqr.gain[0],1,1e-10,"double-integrator gain k1");close(lqr.gain[1],s3,1e-10,"double-integrator gain k2");
    close(lqr.riccati[0],s3,1e-10,"P11");close(lqr.riccati[1],1,1e-10,"P12");close(lqr.riccati[3],s3,1e-10,"P22");
    const auto system=rpe();const auto x=system->controlledInitialState();
    require(x==std::vector<double>({1.,0.}),"controlled equilibrium IC");
    const auto a=system->jacobian(0,x.data()),b=system->inputMatrix(x.data());
    const double nu=system->configuration().dimensionlessGroups.at("viscosity");
    close(a[2],-4*pi*pi,1e-9,"linearized stiffness");close(a[3],-nu,1e-12,"linearized damping");
    const auto klqr=dynsyslab::numerics::solveContinuousLqr(a,b,{10,0,0,1},{1},2,1);
    require(klqr.relativeResidual<1e-10,"Riccati residual");
    for(const auto& pole:klqr.closedLoopPoles)require(pole.real()<0,"stabilizing LQR gain");
}
}

int main(){
    try{checkpointConfiguration();canonicalTrajectory();linearRingdownAndHdmd();riccati();
        std::cout<<"dynsyslab_rpe_hdmd_reference_tests: PASS\n";return 0;}
    catch(const std::exception& e){std::cerr<<"dynsyslab_rpe_hdmd_reference_tests: FAIL: "<<e.what()<<'\n';return 1;}
}
