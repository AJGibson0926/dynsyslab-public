#include "SystemSupport.hpp"
#include <cmath>
#include <stdexcept>

// Canonical Hamilton equations and the Euler Lie-Poisson equations.
// Sources and analytic domains are also attached to AnalyticIntegrable.cpp.
// Arnold, Mathematical Methods of Classical Mechanics, chapters 8-10;
// https://doi.org/10.1007/978-1-4757-2063-1
namespace dynsyslab::core::detail {
namespace {
using E=AnalyticExpression;
double positive(SystemConfiguration& c,const char* key,double fallback) {
    const double v=resolvedParameter(c,key,fallback);
    if(!std::isfinite(v)||v<=0)throw std::invalid_argument(std::string(key)+" must be positive and finite");
    return v;
}
double nonnegative(SystemConfiguration& c,const char* key,double fallback=0) {
    const double v=resolvedParameter(c,key,fallback);
    if(!std::isfinite(v)||v<0)throw std::invalid_argument(std::string(key)+" must be nonnegative and finite");
    return v;
}
void observable(std::vector<Observable>& o,const std::string& key,const E& e) {
    o.push_back(Observable::fromExpressions({key,key,"dimensionless","dimensionless",{},false},{e}));
}
std::unique_ptr<ODESystem> build(SystemConfiguration c,std::vector<double> initial,E H,
    std::vector<E> f,std::vector<E> B,std::vector<Observable> obs={},InvalidFunction invalid={}) {
    const auto d=f.size();std::vector<E> jac;
    for(const auto& v:f)for(std::size_t j=0;j<d;++j)jac.push_back(v.derivative(j));
    auto descriptor=descriptorFor(c.identifier);
    const auto originalCount=obs.size();
    for(std::size_t j=0;j<originalCount;++j) {
        const auto name=obs[j].descriptor().name;
        const auto expression=obs[j].expressions().at(0);
        observable(obs,name+"_dot",expression.flowDerivative(f));
        descriptor.observableTimeDerivatives[name]={name+"_dot",true};
    }
    observable(obs,"total_energy",H);observable(obs,"energy_dot",H.flowDerivative(f));
    observable(obs,"input_power_gain",H.flowDerivative(B));
    for(std::size_t j=0;j<d;++j){const auto name=descriptor.stateNames[j]+"_dot";
        observable(obs,name,f[j]);descriptor.observableTimeDerivatives[descriptor.stateNames[j]]={name,true};}
    descriptor.observableTimeDerivatives["total_energy"]={"energy_dot",true};
    if(!invalid)invalid=[d](double,const double*x){return genericInvalid(x,d);};
    return std::make_unique<FunctionalSystem>(descriptor,c,selectedInitial(c,initial),
        [f,B,d](double,const double*x,double u,double*y){for(std::size_t j=0;j<d;++j)y[j]=f[j].evaluate(x)+u*B[j].evaluate(x);},
        [jac](double,const double*x){std::vector<double> v;for(const auto&e:jac)v.push_back(e.evaluate(x));return v;},
        [B](const double*x){std::vector<double> v;for(const auto&e:B)v.push_back(e.evaluate(x));return v;},
        invalid,std::vector<double>{},std::move(obs),std::map<std::string,std::string>{},
        std::map<std::string,std::vector<std::string>>{{"state",descriptor.stateNames}},f,B);
}
}
std::unique_ptr<ODESystem> makeIntegrableSystem(const SystemConfiguration& configuration) {
    auto c=configuration;const auto& id=c.identifier;
    const bool single=id=="morse_oscillator"||id=="isotonic_oscillator"||id=="poschl_teller"
        ||id=="cosh_oscillator"||id=="radial_twist";
    if(!single&&id!="kepler_planar"&&id!="euler_top"&&id!="coupled_duffing_hamiltonian")return {};
    if(c.preset.empty())c.preset="regular_bounded";
    if(c.preset!="regular_bounded")throw std::invalid_argument("Unknown integrable-system preset");
    const double damping=nonnegative(c,"damping"),gain=resolvedParameter(c,"input_gain",1);
    if(!std::isfinite(gain))throw std::invalid_argument("input_gain must be finite");
    if(single){
        const auto q=E::state(2,0),p=E::state(2,1),zero=E::constant(0,2);
        E H=zero;std::vector<double> initial{.5,0};InvalidFunction invalid;
        if(id=="morse_oscillator"){
            const double D=positive(c,"depth",1),a=positive(c,"inverse_length",1),m=positive(c,"mass",1);
            H=p.power(2)/(2*m)+D*(1-(-a*q).unary("exp")).power(2);
        }else if(id=="isotonic_oscillator"){
            const double w=positive(c,"omega",1),g=positive(c,"inverse_square_strength",1);
            H=p.power(2)/2+w*w*q.power(2)/2+g/(2*q.power(2));initial={1.4,0};
            invalid=[](double,const double*x){return x[0]<=0||genericInvalid(x,2);};
        }else if(id=="poschl_teller"){
            const double D=positive(c,"depth",1),a=positive(c,"inverse_length",1);
            const auto cosh=((a*q).unary("exp")+(-a*q).unary("exp"))/2;
            H=p.power(2)/2-D/cosh.power(2);
        }else if(id=="cosh_oscillator"){
            const double a=positive(c,"strength",1);
            H=p.power(2)/2+a*(q.unary("exp")+(-q).unary("exp"));
        }else{
            const double w=positive(c,"omega",1),alpha=nonnegative(c,"twist",.5);
            const auto I=(q.power(2)+p.power(2))/2;H=w*I+alpha*I.power(2)/2;initial={1,0};
        }
        std::vector<Observable> obs;observable(obs,"position_squared",q.power(2));
        observable(obs,"position_momentum",q*p);
        if(id=="radial_twist") {
            const auto I=(q.power(2)+p.power(2))/2,R=(2*I).power(.5);
            observable(obs,"action",I);
            observable(obs,"hamiltonian_frequency",c.parameters.at("omega")+c.parameters.at("twist")*I);
            // Canonical angle theta=atan2(-p,q); phase coordinates are defined
            // only away from the equilibrium. Their derivatives use full f.
            observable(obs,"phase_cosine",q/R);observable(obs,"phase_sine",-p/R);
        }
        if(id=="poschl_teller") {
            const auto aq=c.parameters.at("inverse_length")*q;
            observable(obs,"hyperbolic_position",(aq.unary("exp")-(-aq).unary("exp"))/2);
        }
        if(id=="cosh_oscillator")
            observable(obs,"half_hyperbolic_position",((q/2).unary("exp")-(-q/2).unary("exp"))/2);
        return build(c,initial,H,{H.derivative(1),-H.derivative(0)-damping*p},{zero,E::constant(gain,2)},obs,invalid);
    }
    if(id=="kepler_planar"){
        const double mu=positive(c,"mu",1);
        const auto q1=E::state(4,0),q2=E::state(4,1),p1=E::state(4,2),p2=E::state(4,3),zero=E::constant(0,4);
        const auto R=(q1.power(2)+q2.power(2)).power(.5),H=(p1.power(2)+p2.power(2))/2-mu/R;
        std::vector<Observable> obs;observable(obs,"angular_momentum",q1*p2-q2*p1);observable(obs,"radius",R);
        observable(obs,"eccentricity_x",(p2*(q1*p2-q2*p1))/mu-q1/R);
        observable(obs,"eccentricity_y",(-p1*(q1*p2-q2*p1))/mu-q2/R);
        return build(c,{1,0,0,.8*std::sqrt(mu)},H,{p1,p2,-H.derivative(0)-damping*p1,-H.derivative(1)-damping*p2},
            {zero,zero,E::constant(gain,4),zero},obs,
            [](double,const double*x){return x[0]*x[0]+x[1]*x[1]<=0||genericInvalid(x,4);});
    }
    if(id=="euler_top"){
        const double I1=positive(c,"inertia_1",1),I2=positive(c,"inertia_2",2),I3=positive(c,"inertia_3",3);
        if(!(I1<I2&&I2<I3))throw std::invalid_argument("Euler top requires inertia_1 < inertia_2 < inertia_3");
        const auto a=E::state(3,0),b=E::state(3,1),z=E::state(3,2),zero=E::constant(0,3);
        const auto H=a.power(2)/(2*I1)+b.power(2)/(2*I2)+z.power(2)/(2*I3);
        std::vector<Observable> obs;observable(obs,"casimir",a.power(2)+b.power(2)+z.power(2));
        return build(c,{.2,.3,1},H,{(1/I3-1/I2)*b*z-damping*a,(1/I1-1/I3)*z*a-damping*b,
            (1/I2-1/I1)*a*b-damping*z},{zero,zero,E::constant(gain,3)},obs);
    }
    const double w1=positive(c,"omega_1",1),w2=positive(c,"omega_2",std::sqrt(2.));
    const double b1=nonnegative(c,"beta_1",1),b2=nonnegative(c,"beta_2",1),eps=nonnegative(c,"coupling",0);
    const auto q1=E::state(4,0),q2=E::state(4,1),p1=E::state(4,2),p2=E::state(4,3),zero=E::constant(0,4);
    const auto H1=p1.power(2)/2+w1*w1*q1.power(2)/2+b1*q1.power(4)/4;
    const auto H2=p2.power(2)/2+w2*w2*q2.power(2)/2+b2*q2.power(4)/4;
    const auto H=H1+H2+eps*q1.power(2)*q2.power(2);
    std::vector<Observable> obs;observable(obs,"component_energy_1",H1);observable(obs,"component_energy_2",H2);
    observable(obs,"coupling_energy",eps*q1.power(2)*q2.power(2));
    return build(c,{.7,.4,0,0},H,{p1,p2,-H.derivative(0)-damping*p1,-H.derivative(1)-damping*p2},
        {zero,zero,E::constant(gain,4),zero},obs);
}
} // namespace dynsyslab::core::detail
