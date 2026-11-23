#include "SystemSupport.hpp"
#include <cmath>
#include <stdexcept>

// Classical model definitions, implemented independently of the class pack.
// Clover: A. Kirichuka (2012), section 4.6, equations (23)-(24), p.97:
// https://www.lumii.lv/uploads/sadirb./2012/LUMII_2012_A_Kirichuka.pdf
// Cubic-quintic multiwell context: https://doi.org/10.1006/jsvi.1995.0448
// No unpublished spectral case-study results are encoded in this file.
namespace dynsyslab::core::detail {
namespace {
void finite(double v) {
    if (!std::isfinite(v)) throw std::invalid_argument("Multiwell parameters must be finite");
}
using Expr = AnalyticExpression;
void add(std::vector<Observable>& out,const char* name,const char* label,const Expr& e) {
    out.push_back(Observable::fromExpressions(
        {name,label,"dimensionless","dimensionless",{},false},{e}));
}
std::unique_ptr<ODESystem> chain(SystemConfiguration c) {
    double a=.5;
    std::vector<double> initial{.2,0};
    if(c.preset.empty()) c.preset="conservative_chain";
    if(c.preset=="left_libration") initial={-1.05,.02};
    else if(c.preset=="right_libration") initial={1.05,.02};
    else if(c.preset=="exterior_orbit") initial={0,.25};
    else if(c.preset=="equal_depth_wells") a=1/std::sqrt(3.);
    else if(c.preset!="conservative_chain") throw std::invalid_argument("Unknown triple-well preset");
    a=resolvedParameter(c,"inner_saddle",a);
    const double damping=resolvedParameter(c,"damping",0);
    const double gain=resolvedParameter(c,"input_gain",1);
    finite(a);finite(damping);finite(gain);
    if(a<=0 || a>=1 || damping<0)
        throw std::invalid_argument("Triple well requires 0<inner_saddle<1 and damping>=0");
    const double b=a*a;
    const auto q=Expr::state(2,0),v=Expr::state(2,1);
    const auto V=q.power(6)/6-(1+b)*q.power(4)/4+b*q.power(2)/2;
    const auto force=-q.power(5)+(1+b)*q.power(3)-b*q-damping*v;
    std::vector<Observable> obs;
    add(obs,"potential_energy","Potential energy",V);
    add(obs,"kinetic_energy","Kinetic energy",v.power(2)/2);
    add(obs,"total_energy","Mechanical energy",V+v.power(2)/2);
    add(obs,"energy_dot","Energy rate (zero input)",-damping*v.power(2));
    add(obs,"input_power_gain","Energy rate per unit input",gain*v);
    add(obs,"force","Acceleration (zero input)",force);
    add(obs,"position_squared","Squared position",q.power(2));
    add(obs,"position_squared_dot","Squared-position rate",2*q*v);
    auto descriptor=descriptorFor("triple_well_duffing");
    return std::make_unique<FunctionalSystem>(descriptor,c,selectedInitial(c,initial),
        [=](double,const double* x,double u,double* f) {
            f[0]=x[1]; f[1]=-x[0]*(x[0]*x[0]-b)*(x[0]*x[0]-1)-damping*x[1]+gain*u;
        },
        [=](double,const double* x) {return std::vector<double>{0,1,
            -5*std::pow(x[0],4)+3*(1+b)*x[0]*x[0]-b,-damping};},
        [=](const double*) {return std::vector<double>{0,gain};},
        [](double,const double* x){return genericInvalid(x,2);},std::vector<double>{0,0},
        std::move(obs),std::map<std::string,std::string>{{"x_1","position"},{"x_2","velocity"}},
        std::map<std::string,std::vector<std::string>>{{"state",{"position","velocity"}},
            {"energy_components",{"potential_energy","kinetic_energy"}}},
        std::vector<Expr>{v,force},std::vector<Expr>{Expr::constant(0,2),Expr::constant(gain,2)});
}
std::unique_ptr<ODESystem> clover(SystemConfiguration c) {
    if(c.preset.empty()) c.preset="right_leaf";
    std::vector<double> initial{.70,.02};
    if(c.preset=="upper_left_leaf" || c.preset=="lower_left_leaf") {
        const double angle=(c.preset=="upper_left_leaf"?2:4)*std::acos(-1.)/3;
        initial={.70*std::cos(angle)-.02*std::sin(angle),.70*std::sin(angle)+.02*std::cos(angle)};
    } else if(c.preset=="exterior_orbit") initial={1.1,0};
    else if(c.preset=="separatrix") initial={1,0};
    else if(c.preset!="right_leaf") throw std::invalid_argument("Unknown clover preset");
    const double gain=resolvedParameter(c,"input_gain",1);finite(gain);
    const auto q=Expr::state(2,0),p=Expr::state(2,1),r2=q.power(2)+p.power(2);
    const auto f1=4*q.power(2)*p+4*p.power(3)+6*q*p;
    const auto f2=3*q.power(2)-4*q.power(3)-3*p.power(2)-4*q*p.power(2);
    const auto H=r2.power(2)-q.power(3)+3*q*p.power(2);
    std::vector<Observable> obs;
    add(obs,"total_energy","Hamiltonian",H);
    add(obs,"energy_dot","Hamiltonian rate (zero input)",Expr::constant(0,2));
    add(obs,"input_power_gain","Hamiltonian rate per unit input",gain*f1);
    add(obs,"x_1_dot","dx_1/dt",f1);
    add(obs,"x_2_dot","dx_2/dt (zero input)",f2);
    add(obs,"radius_squared","Squared phase-plane radius",r2);
    add(obs,"radius_squared_dot","Squared-radius rate (zero input)",2*q*f1+2*p*f2);
    add(obs,"cubic_real","Re((x_1+i x_2)^3)",q.power(3)-3*q*p.power(2));
    add(obs,"cubic_imag","Im((x_1+i x_2)^3)",3*q.power(2)*p-p.power(3));
    add(obs,"hamiltonian_gradient_1","dH/dx_1",-f2);
    add(obs,"hamiltonian_gradient_2","dH/dx_2",f1);
    return std::make_unique<FunctionalSystem>(descriptorFor("three_leaf_hamiltonian_clover"),c,
        selectedInitial(c,initial),
        [=](double,const double* x,double u,double* f) {
            const double q=x[0],p=x[1];
            f[0]=4*q*q*p+4*p*p*p+6*q*p;
            f[1]=3*q*q-4*q*q*q-3*p*p-4*q*p*p+gain*u;
        },
        [](double,const double* x) {const double q=x[0],p=x[1];return std::vector<double>{
            8*q*p+6*p,4*q*q+12*p*p+6*q,6*q-12*q*q-4*p*p,-6*p-8*q*p};},
        [=](const double*) {return std::vector<double>{0,gain};},
        [](double,const double* x){return genericInvalid(x,2);},std::vector<double>{.75,0},
        std::move(obs),std::map<std::string,std::string>{{"x","x_1"},{"y","x_2"}},
        std::map<std::string,std::vector<std::string>>{{"state",{"x_1","x_2"}},
            {"state_and_rate",{"x_1","x_2","x_1_dot","x_2_dot"}},
            {"cubic_coordinates",{"cubic_real","cubic_imag"}}},
        std::vector<Expr>{f1,f2},std::vector<Expr>{Expr::constant(0,2),Expr::constant(gain,2)});
}
}
std::unique_ptr<ODESystem> makeMultiwellSystem(const SystemConfiguration& c) {
    if(c.identifier=="triple_well_duffing")return chain(c);
    if(c.identifier=="three_leaf_hamiltonian_clover")return clover(c);
    return {};
}
}
