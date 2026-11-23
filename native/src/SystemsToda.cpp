#include "SystemSupport.hpp"
#include <cmath>
#include <stdexcept>

// Standard finite Toda chain. Canonical coordinates (q_1,...,q_N,p_1,...,p_N).
// Teschl, Jacobi Operators and Completely Integrable Nonlinear Lattices:
// https://www.mat.univie.ac.at/~gerald/ftp/book-jac/
namespace dynsyslab::core::detail {
std::unique_ptr<ODESystem> makeTodaSystem(const SystemConfiguration& requested) {
    if(requested.identifier!="toda_lattice")return {};
    auto c=requested;if(c.preset.empty())c.preset="periodic";
    if(c.preset!="periodic"&&c.preset!="open")throw std::invalid_argument("Toda preset must be periodic or open");
    const bool periodic=c.preset=="periodic";
    const auto n=resolvedIntegerParameter(c,"sites",3),d=2*n;
    const double alpha=resolvedParameter(c,"alpha",1),damping=resolvedParameter(c,"damping",0);
    const double gain=resolvedParameter(c,"input_gain",1);
    if(n<2||n>32||!std::isfinite(alpha)||alpha<=0||!std::isfinite(damping)||damping<0||!std::isfinite(gain))
        throw std::invalid_argument("Toda requires 2..32 sites, positive alpha, nonnegative damping, finite gain");
    using E=AnalyticExpression;const auto zero=E::constant(0,d);
    std::vector<E> q,p,bonds,f(d,zero),B(d,zero),jac;
    auto H=zero,P=zero,Q=zero,p2=zero;
    auto descriptor=descriptorFor("toda_lattice",d);
    descriptor.stateNames.clear();descriptor.stateLabels.clear();
    std::vector<Observable> observables;
    const auto add=[&](const std::string& name,const E& expression) {
        observables.push_back(Observable::fromExpressions(
            {name,name,"dimensionless","dimensionless",{},false},{expression}));
    };
    for(std::size_t j=0;j<n;++j) {
        q.push_back(E::state(d,j));p.push_back(E::state(d,n+j));
        H=H+p.back().power(2)/2;P=P+p.back();Q=Q+q.back()/double(n);p2=p2+p.back().power(2);
        descriptor.stateNames.push_back("q_"+std::to_string(j+1));
    }
    for(std::size_t j=0;j<n;++j)descriptor.stateNames.push_back("p_"+std::to_string(j+1));
    descriptor.stateLabels=descriptor.stateNames;
    for(std::size_t j=0;j<(periodic?n:n-1);++j) {
        const auto difference=q[j]-q[(j+1)%n];
        bonds.push_back(alpha*alpha*difference.unary("exp"));H=H+bonds.back();
        add("bond_"+std::to_string(j+1),bonds.back());
        add("relative_q_"+std::to_string(j+1),difference);
        add("flaschka_a_"+std::to_string(j+1),(alpha/2)*(difference/2).unary("exp"));
    }
    B[n]=E::constant(gain,d);
    for(std::size_t j=0;j<n;++j) {
        f[j]=p[j];f[n+j]=-damping*p[j];
        if(j>0)f[n+j]=f[n+j]+bonds[j-1];
        else if(periodic)f[n+j]=f[n+j]+bonds[n-1];
        if(j<n-1||periodic)f[n+j]=f[n+j]-bonds[j];
        add("flaschka_b_"+std::to_string(j+1),-p[j]/2);
    }
    for(std::size_t j=0;j<d;++j) {
        add(descriptor.stateNames[j]+"_dot",f[j]);
        descriptor.observableTimeDerivatives[descriptor.stateNames[j]]={descriptor.stateNames[j]+"_dot",j>=n};
        for(std::size_t k=0;k<d;++k)jac.push_back(f[j].derivative(k));
    }
    add("total_energy",H);add("energy_dot",-damping*p2);add("input_power_gain",gain*p[0]);
    add("total_momentum",P);add("momentum_dot",-damping*P);
    add("center_of_mass",Q);add("center_velocity",P/double(n));
    add("relative_energy",H-P.power(2)/(2*double(n)));
    descriptor.observableTimeDerivatives["total_energy"]={"energy_dot",true};
    descriptor.observableTimeDerivatives["total_momentum"]={"momentum_dot",true};
    descriptor.observableTimeDerivatives["center_of_mass"]={"center_velocity",false};
    std::vector<double> initial(d,0);initial[0]=.3;initial[1]=-.3;
    return std::make_unique<FunctionalSystem>(descriptor,c,selectedInitial(c,initial),
        [f,B,d](double,const double*x,double u,double*y){for(std::size_t j=0;j<d;++j)y[j]=f[j].evaluate(x)+u*B[j].evaluate(x);},
        [jac](double,const double*x){std::vector<double> y;for(const auto&e:jac)y.push_back(e.evaluate(x));return y;},
        [B](const double*x){std::vector<double> y;for(const auto&e:B)y.push_back(e.evaluate(x));return y;},
        [d,n](double,const double*x){
            if(genericInvalid(x,d))return true;
            for(std::size_t j=0;j<n;++j)if(std::abs(x[j]-x[(j+1)%n])>600)return true;
            return false;},
        std::vector<double>{},observables,std::map<std::string,std::string>{},
        std::map<std::string,std::vector<std::string>>{{"state",descriptor.stateNames}},f,B);
}
}
