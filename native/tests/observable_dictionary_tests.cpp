#include <dynsyslab/core/BasisFunctions.hpp>
#include <dynsyslab/core/KoopmanSpectralAnalysis.hpp>
#include <dynsyslab/core/System.hpp>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
namespace {
namespace c=dynsyslab::core;
void check(bool value,const char* message){if(!value)throw std::runtime_error(message);}
void close(double a,double b,double tolerance=1e-11) {
    check(std::isfinite(a)&&std::isfinite(b)&&std::abs(a-b)<=tolerance*std::max({1.,std::abs(a),std::abs(b)}),"Value mismatch");
}
template<class F>void rejects(F f){bool failed=false;try{f();}catch(const std::exception&){failed=true;}check(failed,"Expected rejection");}
void observables() {
    for(const std::string system:{"duffing","nonlinear_pendulum"}) {
        c::SystemConfiguration config;config.identifier=system;
        if(system=="nonlinear_pendulum")config.parameters={{"length_m",2},{"gravity_m_s2",8}};
        const auto sys=c::createSystem(config);const auto catalog=sys->observableCatalog();
        double x[]={.7,-1.2};const double v=system=="duffing"?x[1]:.5*x[1];
        const auto g=catalog.at("phase_sech").evaluate(x)[0];
        close(g.real(),std::cos(x[0])/std::cosh(v));close(g.imag(),std::sin(x[0])/std::cosh(v));
        if(system=="duffing")check(std::abs(g-catalog.at("random_observation").evaluate(x)[0])<1e-14,"Legacy observable changed");
        for(const auto& name:catalog.names()) {
            if(name.find("phase")!=0&&name!="gaussian_sum")continue;
            const auto& obs=catalog.at(name);const auto analytic=obs.jacobian(x);
            const auto numeric=obs.jacobian(x,c::ObservableJacobianBackend::finiteDifference);
            for(std::size_t k=0;k<analytic.size();++k)check(std::abs(analytic[k]-numeric[k])<2e-7,"Analytic complex/component Jacobian");
            if(obs.expressions().empty())continue;
            const auto restored=c::Observable::fromExpressionDefinition(c::json::parse(c::json::dump(obs.expressionDefinition())));
            check(restored.evaluate(x)==obs.evaluate(x),"Expression roundtrip");
            const auto f=sys->driftExpressions();const auto rate=obs.flowDerivative(f).evaluate(x);
            for(std::size_t j=0;j<rate.size();++j)
                close(rate[j].real(),(analytic[2*j]*f[0].evaluate(x)+analytic[2*j+1]*f[1].evaluate(x)).real());
            const auto hessian=obs.partialDerivative(0).partialDerivative(1).evaluate(x);
            const auto transposed=obs.partialDerivative(1).partialDerivative(0).evaluate(x);
            for(std::size_t j=0;j<hessian.size();++j)close(hessian[j].real(),transposed[j].real());
        }
        double large[]={.5,1e4};
        for(auto z:catalog.at("phase_sech").evaluate(large))check(std::isfinite(z.real())&&std::isfinite(z.imag()),"Stable sech");
    }
}
void bases() {
    const double x=.37,pi=3.14159265358979323846;
    close(c::basis::chebyshev(4,x),8*std::pow(x,4)-8*x*x+1);
    close(c::basis::laguerre(3,x),1-3*x+1.5*x*x-x*x*x/6);
    const double h0=std::pow(pi,-.25)*std::exp(-x*x/2);
    close(c::basis::hermiteFunction(2,x),(2*x*x-1)*h0/std::sqrt(2.));
    check(c::basis::hermiteFunction(512,40)>0,"Hermite Gaussian underflow must not erase high modes");
    close(c::basis::radialKernel("matern",.7,1,1.5),1.7*std::exp(-.7));
    close(c::basis::radialKernel("matern",0,1,2,2),2);
    close(c::basis::radialKernel("matern",1e-100,1,8),1);
    rejects([]{c::basis::radialKernel("matern",1,0);});
    // Quadrature checks normalization and orthogonality, not just recurrence identities.
    double inner=0,norm=0;const double step=.005;
    for(int j=0;j<=4000;++j) {
        const double t=-10+j*step,w=(j==0||j==4000)?.5:1;
        inner+=w*step*c::basis::hermiteFunction(2,t)*c::basis::hermiteFunction(4,t);
        norm+=w*step*std::pow(c::basis::hermiteFunction(4,t),2);
    }
    close(inner,0);close(norm,1);
}
void dictionaries() {
    const std::vector<double> points{0,0,1,1,-1,-1,.2,.3,2,-2};
    c::SpectralDictionaryOptions o;o.kind="chebyshev";o.degree=9;o.indexSet="tensor_product";o.scale={2.5,2.5};
    auto d=c::buildSpectralDictionary(points,2,o);check(d.size()==100,"Chebyshev reference size");
    auto restored=c::spectralDictionaryFromJson(c::json::parse(c::json::dump(c::spectralDictionaryJson(d))));
    check(restored.evaluate(points)==d.evaluate(points),"Frozen basis roundtrip");
    o={};o.kind="laguerre";o.degree=12;d=c::buildSpectralDictionary(points,2,o);check(d.size()==91,"Laguerre reference size");
    o={};o.kind="fourier_hermite";o.angularOrder=3;o.hermiteOrder=4;
    d=c::buildSpectralDictionary(points,2,o);check(d.size()==35,"Fourier-Hermite size");
    const auto a=d.evaluate({.7,-.3}),b=d.evaluate({.7+2*3.14159265358979323846,-.3});
    for(std::size_t j=0;j<a.size();++j)close(a[j],b[j]);
    o.normalize=true;rejects([&]{c::buildSpectralDictionary(points,2,o);});
    o={};o.kind="polynomial";o.degree=32;o.maximumFeatures=10;
    rejects([&]{c::buildSpectralDictionary(points,2,o);});
    o={};o.kind="laplacian";o.centers=2;o.centerSelection="kmeans";
    d=c::buildSpectralDictionary({-2,-2,0,0,8,8,10,10},2,o);
    check(d.nodes==std::vector<double>({-1,-1,9,9}),"Deterministic kmeans centers");
    for(const std::string kind:{"linear","polynomial","chebyshev","laguerre","hermite","gaussian","laplacian","lorentzian","matern","partition"}) {
        o={};o.kind=kind;o.centers=4;d=c::buildSpectralDictionary(points,2,o);
        restored=c::spectralDictionaryFromJson(c::json::parse(c::json::dump(c::spectralDictionaryJson(d))));
        check(restored.evaluate(points)==d.evaluate(points),"Dictionary roundtrip");
    }
}
}
int main() {
    try{observables();bases();dictionaries();std::cout<<"Observable/dictionary tests passed\n";return 0;}
    catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 1;}
}
