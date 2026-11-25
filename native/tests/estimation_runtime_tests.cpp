#include <dynsyslab/core/EstimationRuntime.hpp>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
namespace c=dynsyslab::core;
void check(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
void close(double a,double b,double tol=1e-11){check(std::abs(a-b)<tol,"Observer value differs from analytic reference");}
template<class F> void rejects(F f) {
    try{f();}catch(const std::logic_error&){return;}
    throw std::runtime_error("Invalid observer input was accepted");
}
c::KoopmanModel scalar() {
    c::FunctionDomain q;q.dimension=1;
    c::EigenfunctionFitResult bank{c::FunctionExpansion({1,{1},{0},{1}}, {1.}, q), c::Spectrum::fromContinuousRates({-1.}), {}};
    c::KoopmanModel model(bank,1);
    model.setControlMap("fit",c::ConstantControlMap{{2.},1,1,"test fixture"});
    model.setReconstruction("y",{{1.},1,1,0,0,0});
    model.setReconstruction("x",{{3.},1,1,0,0,0});
    return model;
}
void heldScalar() {
    auto model=scalar();
    const auto design=c::designKoopmanLQEFromMeasurements(model,"y",{1},{1},"x");
    const double L=design.solution.gain[0],a=-1-L,h=.7;
    c::KoopmanLQEObserver observer(model,design,"fit",{.4},{10},{20});
    close(observer.observe(5,{13},{0}).output[0],21.2);
    const auto value=observer.observe(5+h,{100},{.2});
    const double expected=std::exp(a*h)*.4+std::expm1(a*h)/a*(L*3+2*.2);
    close(value.lifted[0],expected);close(value.output[0],20+3*expected);
    close(value.measurement[0],10+expected);close(value.innovation[0],90-expected);
    // Current measurement has no retroactive influence; caller supplies actual u.
    rejects([&]{observer.observe(5+h,{0},{0});});
    rejects([&]{observer.observe(5+h+1,{0},{});});
    rejects([&]{observer.observe(5+h+1,{std::numeric_limits<double>::quiet_NaN()},{0});});
    close(observer.estimate().lifted[0],expected);
    model.setControlMap("fit",c::ConstantControlMap{{99.},1,1,"later mutation"});
    const auto next=observer.observe(5+2*h,{10},{.2});
    close(next.lifted[0],std::exp(a*h)*expected+std::expm1(a*h)/a*(L*90+.4));
    rejects([&]{c::KoopmanLQEObserver bad(model,design,"exact",{0});});
    auto corrupt=design;corrupt.measurementMatrix[0]=2;
    rejects([&]{c::KoopmanLQEObserver bad(model,corrupt,"fit",{0});});
    model.setDiscreteInputMap("discrete",{{1.},1,1,.1,"sampled"});
    rejects([&]{c::KoopmanLQEObserver bad(model,design,"discrete",{0});});
}
void autonomousParity() {
    const auto model=scalar();
    const auto design=c::designKoopmanLQEFromMeasurements(model,"y",{1},{.2},"x");
    const std::vector<double> t{0,.1,.2,.3,.5},y{1,.8,.6,.2,.1};
    const auto batch=c::runKoopmanLQEFromLifted(model,design,t,y,{.3});
    c::KoopmanLQEObserver online(model,design,"fit",{.3});
    for(std::size_t k=0;k<t.size();++k) {
        const auto sample=online.observe(t[k],{y[k]},{0});
        close(sample.lifted[0],batch.liftedEstimate[k]);
        close(sample.output[0],batch.stateEstimate[k]);
    }
}
void measurementInputTiming() {
    const auto model=scalar();
    const auto design=c::designKoopmanLQEFromMeasurements(model,"y",{1},{1},"x");
    const double L=design.solution.gain[0],h=.2,a=-1-L,D=3;
    c::KoopmanLQEObserver observer(model,design,"fit",{.4});
    observer.setMeasurementFeedthrough({D});
    // y(t0) was sampled before the controller changed input from .2 to .8.
    observer.observe(0,{.7+D*.2},{.2});
    const auto next=observer.observe(h,{1.1+D*.8},{.8});
    const double expected=std::exp(a*h)*.4+std::expm1(a*h)/a*(L*.7+2*.8);
    close(next.lifted[0],expected);close(next.measurement[0],expected+D*.8);
    close(next.innovation[0],1.1-expected);
    c::KoopmanLQEObserver scaled(model,design,"fit",{.4},{10},{20});
    scaled.setMeasurementFeedthrough({D},0);
    double z=.4;
    const double oldY=10+.7+D*.2*(20+3*.7),newY=10+1.1+D*.8*(20+3*1.1);
    scaled.observe(0,{oldY},{.2});
    const double scaledA=-1-L*(1+3*D*.2);
    z=std::exp(scaledA*h)*z+std::expm1(scaledA*h)/scaledA*(L*(oldY-10-D*.2*20)+2*.8);
    const auto first=scaled.observe(h,{newY},{.8});
    close(first.lifted[0],z);close(first.measurement[0],10+z+D*.8*(20+3*z));
    // Actuation is unchanged, but the input associated with the held sensor
    // sample changed. Its observation matrix must invalidate the ZOH cache.
    const double nextA=-1-L*(1+3*D*.8);
    z=std::exp(nextA*h)*z+std::expm1(nextA*h)/nextA*(L*(newY-10-D*.8*20)+2*.8);
    close(scaled.observe(2*h,{newY},{.8}).lifted[0],z);
    rejects([&]{scaled.setMeasurementFeedthrough({D},-2);});
}
void bilinear() {
    auto model=scalar();
    model.setControlMap("bilinear",c::BilinearControlMap{{.5},{2.},{},1,1,"test"});
    const auto design=c::designKoopmanLQEFromMeasurements(model,"y",{1},{1});
    c::KoopmanLQEObserver observer(model,design,"bilinear",{.4});
    observer.observe(0,{.7},{0});
    double z=.4,t=0;
    for(double u:{.2,-.6,1.1}) {
        const double a=-1-design.solution.gain[0]+.5*u,h=.2;
        z=std::exp(a*h)*z+std::expm1(a*h)/a*(design.solution.gain[0]*.7+2*u);
        t+=h;close(observer.observe(t,{.7},{u}).lifted[0],z);
    }
    // Complex/conjugate map: compare to independent small-step complex RK4.
    c::FunctionDomain q;q.dimension=2;
    c::EigenfunctionFitResult bank{c::FunctionExpansion({2,{1,0,0,1},{0,0},{1,1}}, {1.,{0,1}}, q), c::Spectrum::fromContinuousRates({{-.3,2}}), {}};
    c::KoopmanModel oscillator(bank,1);
    oscillator.setReconstruction("y",{{1.},1,1,0,0,0});
    const std::complex<double> N{.2,.1},M{.1,-.2},b{.3,.4};
    oscillator.setControlMap("b",c::BilinearControlMap{{N},{b},{M},1,1,"test"});
    const auto d=c::designKoopmanLQEFromMeasurements(oscillator,"y",{1,0,0,1},{.2});
    c::KoopmanLQEObserver o(oscillator,d,"b",{.4,-.2});
    o.observe(0,{.8},{0});const auto actual=o.observe(.2,{.1},{.7}).lifted;
    std::complex<double> reference{.4,-.2},gain{d.solution.gain[0],d.solution.gain[1]};
    const auto f=[&](auto z){return std::complex<double>(-.3,2)*z
        +.7*(b+N*z+M*std::conj(z))+gain*(.8-z.real());};
    for(int k=0;k<2000;++k) {
        const double h=.0001;const auto a=f(reference),b1=f(reference+h*a/2.);
        const auto c1=f(reference+h*b1/2.),d1=f(reference+h*c1);
        reference+=h*(a+2.*b1+2.*c1+d1)/6.;
    }
    close(actual[0],reference.real());close(actual[1],reference.imag());
    check(oscillator.splitValues(oscillator.valuesFromSplit(actual))==actual,"Split-real round trip failed");
}
}
int main(){try{
    heldScalar();autonomousParity();bilinear();measurementInputTiming();
    std::cout<<"Input-aware online KLQE tests passed\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
