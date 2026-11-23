#include <dynsyslab/core/System.hpp>
#include <dynsyslab/core/Toda.hpp>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace dynsyslab::core;
void check(bool b,const char* m){if(!b)throw std::runtime_error(m);}
void verify(std::size_t n,bool periodic) {
    SystemConfiguration c;c.identifier="toda_lattice";c.preset=periodic?"periodic":"open";
    c.parameters={{"sites",double(n)},{"alpha",1.3}};
    auto s=createSystem(c);std::vector<double> x(2*n),f(2*n);
    for(std::size_t j=0;j<n;++j){x[j]=.3*std::sin(double(j+1));x[n+j]=.2*std::cos(double(j+1));}
    s->rhsWithInput(0,x.data(),0,f.data());
    const auto H=s->observableCatalog().at("total_energy");
    const auto grad=H.jacobian(x.data());double work=0,momentum=0;
    for(std::size_t j=0;j<2*n;++j)work+=grad[j].real()*f[j];
    for(std::size_t j=0;j<n;++j)momentum+=f[n+j];
    check(std::abs(work)<1e-12&&std::abs(momentum)<1e-12,"Toda energy and momentum invariants");
    auto pair=todaLaxPair({x.begin(),x.begin()+n},{x.begin()+n,x.end()},1.3,periodic);
    auto xm=x,xp=x;const double h=1e-6;
    for(std::size_t j=0;j<2*n;++j){xp[j]+=h*f[j];xm[j]-=h*f[j];}
    auto plus=todaLaxPair({xp.begin(),xp.begin()+n},{xp.begin()+n,xp.end()},1.3,periodic);
    auto minus=todaLaxPair({xm.begin(),xm.begin()+n},{xm.begin()+n,xm.end()},1.3,periodic);
    for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j) {
        double commutator=0;for(std::size_t k=0;k<n;++k)
            commutator+=pair.generator[i*n+k]*pair.matrix[k*n+j]-pair.matrix[i*n+k]*pair.generator[k*n+j];
        check(std::abs((plus.matrix[i*n+j]-minus.matrix[i*n+j])/(2*h)-commutator)<1e-8,"Toda Lax commutator");
    }
    auto jac=s->jacobian(0,x.data());
    for(std::size_t j=0;j<2*n;++j) {
        auto a=x,b=x;a[j]+=h;b[j]-=h;std::vector<double> fa(2*n),fb(2*n);
        s->rhsWithInput(0,a.data(),0,fa.data());s->rhsWithInput(0,b.data(),0,fb.data());
        for(std::size_t i=0;i<2*n;++i)check(std::abs(jac[i*2*n+j]-(fa[i]-fb[i])/(2*h))<1e-8,"Toda Jacobian");
    }
    auto family=s->analyticSolutions().at(0);
    if(n!=2){check(!family.instantiate(x,0),"General theta evaluator must not be claimed");return;}
    auto orbit=family.instantiate(x,1.7);check(bool(orbit),"Two-site analytic solution");
    for(double t:{-8.,-.2,1.7,2.3,6.1,17.}) {
        auto y=orbit.evaluate(t),a=orbit.evaluate(t+h),b=orbit.evaluate(t-h);
        s->rhsWithInput(t,y.data(),0,f.data());
        for(std::size_t j=0;j<4;++j)check(std::abs((a[j]-b[j])/(2*h)-f[j])<1e-7,"Toda analytic physical-time ODE");
        check(std::abs(H.evaluate(y.data())[0]-H.evaluate(x.data())[0])<1e-10,"Toda analytic energy");
    }
    check(!orbit.torus&&orbit.period==0,"Moving center of mass is not a periodic state");
    x[3]=-x[2];orbit=family.instantiate(x,0);
    if(periodic){check(bool(orbit.torus)&&orbit.period>0,"Zero-momentum periodic Toda torus");
        auto a=orbit.evaluate(orbit.period);
        for(std::size_t j=0;j<4;++j)check(std::abs(a[j]-x[j])<1e-9,"Toda period");}
    else check(!orbit.torus&&orbit.period==0,"Open Toda scattering is not a torus");
}
int main(){try {
    for(bool periodic:{false,true})for(std::size_t n:{2,3,4,8})verify(n,periodic);
    std::cout<<"Toda system, Lax and analytic tests passed\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
