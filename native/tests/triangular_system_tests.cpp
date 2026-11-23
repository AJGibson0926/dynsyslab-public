#include <dynsyslab/core/KnownEigenfunctions.hpp>
#include <dynsyslab/core/Integrator.hpp>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace dynsyslab::core;
namespace {
void check(bool v,const char* label) { if(!v) throw std::runtime_error(label); }
void near(double x,double y,double tol=1e-12) {
    check(std::isfinite(x)&&std::isfinite(y)&&std::abs(x-y)<=tol*(1+std::abs(y)),
          "Triangular numeric identity failed");
}
std::vector<double> product(const std::vector<double>& a,const std::vector<double>& b) {
    std::vector<double> c(16);
    for(size_t i=0;i<4;++i) for(size_t j=0;j<4;++j) for(size_t k=0;k<4;++k)
        c[4*i+j]+=a[4*i+k]*b[4*k+j];
    return c;
}
void verify(SystemConfiguration config) {
    auto sys=createSystem(config);
    const auto& p=sys->configuration().parameters;
    auto exact=knownPolynomialEigenfunctions(*sys);
    check(exact.eigenfunctions.has_value(),"Exact bank missing");
    const auto& bank=*exact.eigenfunctions;
    const auto rates=bank.spectrum().continuousRates();
    check(bank.eigenfunctionCount()==4 && bank.coordinates().names[0]=="x_1","Bank contract");
    check(sys->descriptor().stateNames==std::vector<std::string>{"x_1","x_2"},"State names");
    check(sys->observableCatalog().contains("eta"),"Analytic observable missing");
    for(double x1:{-2.,-.1,0.,.7,3.}) for(double x2:{-1.,0.,2.}) {
        const double x[]{x1,x2}; double f[2];
        sys->rhsWithInput(0,x,0,f);
        auto z=bank.evaluate(x);
        for(size_t m=0;m<4;++m) {
            const auto gradient=bank.gradient(m,x);
            near((gradient[0]*f[0]+gradient[1]*f[1]).real(),(rates[m]*z[m]).real());
            for(size_t j=0;j<2;++j) {
                std::complex<double> expected{};
                for(size_t n=0;n<4;++n) expected+=exact.derivativeOperators[16*j+4*m+n]*z[n];
                near(gradient[j].real(),expected.real());
            }
        }
        for(size_t j=0;j<2;++j) {
            std::complex<double> reconstructed{};
            for(size_t n=0;n<4;++n) reconstructed+=exact.stateReadout[4*j+n]*z[n];
            near(reconstructed.real(),x[j]);
        }
        const auto jac=sys->jacobian(0,x);
        const double eps=1e-5;
        for(size_t j=0;j<2;++j) {
            double xp[]{x1,x2},xm[]{x1,x2},fp[2],fm[2];
            xp[j]+=eps; xm[j]-=eps;
            sys->rhsWithInput(0,xp,0,fp); sys->rhsWithInput(0,xm,0,fm);
            for(size_t i=0;i<2;++i) near(jac[2*i+j],(fp[i]-fm[i])/(2*eps),1e-9);
        }
        const auto g=sys->inputMatrix(x);
        double controlled[2]; sys->rhsWithInput(0,x,.37,controlled);
        for(size_t j=0;j<2;++j) near(controlled[j]-f[j],.37*g[j]);
    }
    std::vector<double> d1(exact.derivativeOperators.begin(),exact.derivativeOperators.begin()+16);
    std::vector<double> d2(exact.derivativeOperators.begin()+16,exact.derivativeOperators.end());
    check(product(d1,d2)==product(d2,d1),"Mixed derivatives do not commute");
    check(product(product(d1,d1),d1)==std::vector<double>(16),"Third derivative did not vanish");
    check(product(d2,d2)==std::vector<double>(16),"Second x2 derivative did not vanish");
    auto family=sys->analyticSolutions().at(0);
    const std::vector<double> initial{1.2,-.3};
    auto flow=family.instantiate(initial,.371);
    check(bool(flow),"Exact trajectory missing");
    IntegrationOptions options;
    options.initialTime=.371; options.finalTime=3.371; options.outputStep=.037;
    options.relativeTolerance=1e-11; options.absoluteTolerance=1e-13;
    auto numerical=integrate(*sys,initial,1,options);
    for(size_t n=0;n<numerical.sampleCount;++n) {
        const auto expected=flow.evaluate(numerical.time[n]);
        for(size_t j=0;j<2;++j) near(numerical.state(0,n)[j],expected[j],2e-9);
    }
    const auto json=knownEigenfunctionsJson(exact);
    check(json.at("mode_count").number()==4 && json.at("derivative_operators").array().size()==32,
          "Exact bank JSON lost shape");
    const double a=p.at("a"),b=p.at("b"),c=p.at("c"),k=c==0?0:c/(b-2*a);
    for(double t:{-.5,0.,.4,4.,1000.}) {
        const auto x=flow.evaluate(t+.371);
        near(x[1]+k*x[0]*x[0],std::exp(b*t)*(initial[1]+k*initial[0]*initial[0]),2e-11);
    }
}
}
int main() {
    try {
        SystemConfiguration c; c.identifier="quadratic_triangular";
        verify(c);
        for(const char* preset:{"brunton_2016","kaiser_2021"}) { c.preset=preset; verify(c); }
        c.preset="gradient_closed_reference"; c.parameters={{"input_gain_1",.3},{"input_gain_2",-.7}};
        verify(c);
        c.parameters={{"a",-1},{"b",-2},{"c",1}};
        auto resonant=createSystem(c);
        check(!knownPolynomialEigenfunctions(*resonant).eigenfunctions,"Resonance silently diagonalized");
        auto flow=resonant->analyticSolutions().at(0).instantiate({2,1},.7);
        const auto x=flow.evaluate(1.2);
        near(x[0],2*std::exp(-.5)); near(x[1],3*std::exp(-1.));
        (void)resonant->observableCatalog();
        c.parameters["c"]=0; verify(c); // resonance is harmless when uncoupled
        c.parameters={{"b",-2+1e-10}};
        auto almost=createSystem(c)->analyticSolutions().at(0).instantiate({2,1},0);
        near(almost.evaluate(.5)[1],3*std::exp(-1.),1e-9);
        c.parameters.clear(); c.driver.amplitude=.1;
        check(!createSystem(c)->analyticSolutions().at(0).instantiate({1,0},0),
              "Autonomous formula accepted forcing");
        c.driver.amplitude=0; c.parameters["a"]=std::numeric_limits<double>::quiet_NaN();
        bool rejected=false;
        try {(void)createSystem(c);} catch(const std::invalid_argument&) {rejected=true;}
        check(rejected,"Nonfinite parameter accepted");
        c.parameters.clear(); c.preset="typo"; rejected=false;
        try {(void)createSystem(c);} catch(const std::invalid_argument&) {rejected=true;}
        check(rejected,"Unknown preset silently selected");
        c.preset.clear();
        auto reference=createSystem(c);
        for(double C:{-2.,0.,.7}) {
            const double x1=1.3;
            auto curve=reference->analyticSolutions().at(0).instantiate(
                {x1,x1*x1+C*x1*x1*x1},0);
            for(double t:{0.,.3,2.}) {
                const auto p=curve.evaluate(t);
                near(p[1],p[0]*p[0]+C*p[0]*p[0]*p[0]);
            }
        }
        check(reference->dynamicalKnowledge().papers.size()==2,"Paper provenance missing");
        c.preset="infinite_closure_saddle";
        auto saddle=createSystem(c);
        const auto catalog=saddle->observableCatalog();
        auto seed=catalog.at("invariant_exponential").expressions().at(0);
        auto derivative=seed;
        const double sample[]{.7,-.2};
        for(int n=0;n<=6;++n) {
            near(derivative.evaluate(sample),
                 std::pow(sample[0],3*n)*seed.evaluate(sample),2e-12);
            near(derivative.flowDerivative(saddle->driftExpressions()).evaluate(sample),
                 3*n*derivative.evaluate(sample),2e-11);
            derivative=derivative.derivative(1);
        }
        near(seed.flowDerivative(saddle->driftExpressions()).evaluate(sample),0,1e-13);
        check(bool(knownPolynomialEigenfunctions(*saddle).eigenfunctions),
              "Saddle polynomial bank should still close");
        std::cout<<"triangular system: exact flow, resonance, eigenfunctions, derivatives, "
                    "readout, input and invariant-curve checks passed\n";
        return 0;
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n'; return 1;}
}
