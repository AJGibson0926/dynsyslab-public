#include "AnalyticSystems.hpp"
#include <dynsyslab/core/SpecialFunctions.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>

// Independently derived Hamiltonian solutions, with physical-time phases.
// Morse: Krajňák & Wiggins https://arxiv.org/abs/1905.04059;
// DeMarcus https://doi.org/10.1119/1.11110 (also Fourier coefficients).
// Isotonic: Pinney, Proc. AMS 1 (1950), 681;
// https://doi.org/10.1090/S0002-9939-1950-0038586-4.
// Poeschl-Teller: Blazquez & Negro, J. Phys. A 57 (2024) 195204.
// Elliptic identities: https://dlmf.nist.gov/22.13, 23.6.
namespace dynsyslab::core::detail {
std::vector<AnalyticSolutionFamily> celestialAnalyticSolutions(const SystemConfiguration&);
namespace {
constexpr double pi=3.14159265358979323846;
double p(const SystemConfiguration&c,const char*k){return c.parameters.at(k);}
AnalyticTrajectory no(const char*r){return {{},r};}
AnalyticTrajectory circle(std::function<AnalyticState(double)> map,double phase,double w,double t0){
    AnalyticTrajectory a{[=](double t){return map(phase+w*(t-t0));},{},2*pi/w};
    a.angularFrequencies={w};a.torus=[=](const std::vector<double>&theta){
        if(theta.size()!=1||!std::isfinite(theta[0]))throw std::invalid_argument("Expected one finite torus angle");
        return map(phase+theta[0]);};return a;
}
}
std::vector<AnalyticSolutionFamily> integrableAnalyticSolutions(const SystemConfiguration& c){
    const auto& id=c.identifier;
    if(id=="toda_lattice")return todaAnalyticSolutions(c);
    const bool single=id=="morse_oscillator"||id=="isotonic_oscillator"||id=="poschl_teller"
        ||id=="cosh_oscillator"||id=="radial_twist";
    if(!single&&id!="coupled_duffing_hamiltonian")return celestialAnalyticSolutions(c);
    AnalyticSolutionDescriptor d{id+"_regular",id+" regular invariant set",
        "\\dot\\theta=\\Omega,\\quad x=K(\\theta)",
        "Zero damping and driver; regular bounded orbit; coupled pair requires zero coupling",
        "elementary or Jacobi/Weierstrass evaluation in model time",
        {"https://doi.org/10.1007/978-1-4757-2063-1","https://dlmf.nist.gov/22.13"}};
    if(id=="morse_oscillator")d.references={"https://arxiv.org/abs/1905.04059","https://doi.org/10.1119/1.11110"};
    if(id=="isotonic_oscillator")d.references={"https://doi.org/10.1090/S0002-9939-1950-0038586-4"};
    if(id=="poschl_teller")d.references={"https://uvadoc.uva.es/handle/10324/81442"};
    return {{d,[=](const AnalyticState& x,double t0)->AnalyticTrajectory{
        if(c.driver.amplitude!=0||p(c,"damping")!=0)return no("Regular invariant-torus evaluator requires zero forcing and damping");
        if(id=="coupled_duffing_hamiltonian"){
            if(p(c,"coupling")!=0)return no("Nonzero coupling is Hamiltonian but separable integrability is not established");
            const double w1=p(c,"omega_1"),w2=p(c,"omega_2"),b1=p(c,"beta_1"),b2=p(c,"beta_2");
            const double E1=x[2]*x[2]/2+w1*w1*x[0]*x[0]/2+b1*std::pow(x[0],4)/4;
            const double E2=x[3]*x[3]/2+w2*w2*x[1]*x[1]/2+b2*std::pow(x[1],4)/4;
            if(E1<=0||E2<=0)return no("This torus evaluator requires both component energies positive");
            const special::QuarticOrbit a({2*E1,0,-w1*w1,0,-b1/2},x[0],x[2]);
            const special::QuarticOrbit b({2*E2,0,-w2*w2,0,-b2/2},x[1],x[3]);
            const double o1=2*pi/a.realPeriod(),o2=2*pi/b.realPeriod();
            AnalyticTrajectory result{[=](double t){auto aa=a.evaluate(t-t0),bb=b.evaluate(t-t0);return AnalyticState{aa[0],bb[0],aa[1],bb[1]};},{}};
            result.angularFrequencies={o1,o2};result.torus=[=](const std::vector<double>&theta){
                if(theta.size()!=2||!std::isfinite(theta[0])||!std::isfinite(theta[1]))throw std::invalid_argument("Expected two finite torus angles");
                auto aa=a.evaluate(theta[0]/o1),bb=b.evaluate(theta[1]/o2);return AnalyticState{aa[0],bb[0],aa[1],bb[1]};};return result;
        }
        const double q=x[0],v=x[1];
        if(id=="morse_oscillator"){
            const double D=p(c,"depth"),a=p(c,"inverse_length"),m=p(c,"mass");
            const double E=v*v/(2*m)+D*std::pow(1-std::exp(-a*q),2);
            if(!(E>0&&E<D))return no("Morse periodic annulus requires 0 < energy < depth");
            const double k=std::sqrt(E/D),w=a*std::sqrt(2*(D-E)/m),y=std::exp(a*q);
            const double phase=std::atan2(-a*v*y*(1-E/D)/(m*k*w),(y*(1-E/D)-1)/k);
            return circle([=](double z){const double c=1+k*std::cos(z);
                return AnalyticState{std::log(c/(1-E/D))/a,-m*k*w*std::sin(z)/(a*c)};},phase,w,t0);
        }
        if(id=="isotonic_oscillator"){
            if(q<=0)return no("Isotonic positive component requires q > 0");
            const double w=p(c,"omega"),g=p(c,"inverse_square_strength");
            const double E=(v*v+w*w*q*q+g/(q*q))/2,delta=E*E-g*w*w;
            if(delta<=1e-14*E*E)return no("Equilibrium has no regular phase circle");
            const double A=std::sqrt(delta)/ (w*w),phase=std::atan2(-q*v/w,(q*q-E/(w*w)));
            return circle([=](double z){double u=std::sqrt(E/(w*w)+A*std::cos(z));
                return AnalyticState{u,-w*A*std::sin(z)/u};},phase,2*w,t0);
        }
        if(id=="poschl_teller"){
            const double D=p(c,"depth"),a=p(c,"inverse_length"),E=v*v/2-D/std::pow(std::cosh(a*q),2);
            if(!(E>-D&&E<0))return no("Poeschl-Teller periodic annulus requires -depth < energy < 0");
            const double A=std::sqrt((E+D)/(-E)),w=a*std::sqrt(-2*E);
            const double phase=std::atan2(std::sinh(a*q)/A,a*std::cosh(a*q)*v/(A*w));
            return circle([=](double z){double u=A*std::sin(z);return AnalyticState{std::asinh(u)/a,A*w*std::cos(z)/(a*std::sqrt(1+u*u))};},phase,w,t0);
        }
        if(id=="radial_twist"){
            const double R=std::hypot(q,v),I=R*R/2,w=p(c,"omega")+p(c,"twist")*I;
            if(R==0)return no("Equilibrium has no regular phase circle");
            return circle([=](double z){return AnalyticState{R*std::cos(z),-R*std::sin(z)};},std::atan2(-v,q),w,t0);
        }
        const double a=p(c,"strength"),E=v*v/2+2*a*std::cosh(q),u=std::sinh(q/2),udot=.5*std::cosh(q/2)*v;
        if(E<=2*a)return no("Cosh equilibrium has no regular phase circle");
        const special::QuarticOrbit orbit({(E-2*a)/2,0,(E-6*a)/2,0,-2*a},u,udot);
        const double T=orbit.realPeriod();
        return circle([=](double phase){auto z=orbit.evaluate(phase*T/(2*pi));return AnalyticState{2*std::asinh(z[0]),2*z[1]/std::sqrt(1+z[0]*z[0])};},0,2*pi/T,t0);
    }}};
}
} // namespace
