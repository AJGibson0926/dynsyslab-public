#include <dynsyslab/core/EnsembleBayes.hpp>
#include <dynsyslab/core/PeriodicEstimation.hpp>
#include <dynsyslab/math/PeriodicAtlas.hpp>
#include <dynsyslab/numerics/NonlinearLeastSquares.hpp>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace dynsyslab;
void check(bool x){if(!x)throw std::runtime_error("Periodic estimation test failed");}
void near(double a,double b,double tol=1e-10){check(std::abs(a-b)<tol*(1+std::abs(b)));}
template<class F>void rejects(F f){bool caught=false;try{f();}catch(const std::exception&){caught=true;}check(caught);}
int main(){try{
    core::EnsembleBayesFilter f({{-1},{1}});
    const auto id=[](const auto& x){return x;};
    near(f.effectiveSampleSize(),2);f.correct({1},{1},id);
    near(f.weights()[1],1/(1+std::exp(-2.)));near(f.mean(id)[0],std::tanh(1.));
    f.predict([](const auto& x){return std::vector<double>{x[0]+3};});near(f.mean(id)[0],3+std::tanh(1.));
    const auto before=f.weights();rejects([&]{f.correct({0},{0},id);});check(before==f.weights());
    rejects([]{core::EnsembleBayesFilter bad({{0},{1}},{0,0});});
    math::FourierSeries a{1,{{-1,.5},{1,.5}}},b{2,{{-1,1.},{1,1.}}};
    math::PeriodicAtlas atlas({{1,{a}},{2,{b}}});
    near(atlas.evaluate(1.5,.3)[0].real(),1.5*std::cos(.3));
    near(atlas.predict(1.5,.3,2)[0].real(),1.5*std::cos(3.3));
    const math::PeriodicAtlasTable table(atlas,4096);
    near(table.frequency(1.4823),1.4823);
    near(table.evaluate(1.4823,.371)[0].real(),1.4823*std::cos(.371),1e-6);
    rejects([&]{atlas.evaluate(.9,0);});
    const auto values=math::uniformFourierSamples(a,63);
    for(int j=0;j<63;++j)near(values[j].real(),std::cos(2*std::acos(-1.)*j/63));
    // A genuinely held-out parameter, present in the filter search but absent
    // from training knots; unknown initial phase represented by a uniform prior.
    core::PeriodicFilterOptions o;o.parameters={1,1.5,2};o.phaseCount=32;o.tableSamples=1024;
    o.measurementChannels={0};o.outputChannels={0};o.measurementVariances={.01};
    std::vector<double> time,y;const double phase=2*std::acos(-1.)*7/32;
    for(int j=0;j<80;++j){time.push_back(j*.1);y.push_back(1.5*std::cos(phase+1.5*time.back()));}
    const auto result=core::filterPeriodicAtlas(atlas,time,y,o);
    near(result.parameterMean.back(),1.5,1e-8);near(result.outputs.back(),y.back(),2e-5);
    // Causality: appending observations does not change previous outputs.
    const auto shortRun=core::filterPeriodicAtlas(atlas,{time.begin(),time.begin()+20},{y.begin(),y.begin()+20},o);
    for(int j=0;j<20;++j)near(result.outputs[j],shortRun.outputs[j]);
    std::vector<double> continuous;
    for(double t:time)continuous.push_back(1.4823*std::cos(.371+1.4823*t));
    o.tableSamples=8192;
    const auto identified=core::identifyPeriodicOrbit(atlas,time,continuous,{{1.5,.4}},o);
    near(identified.parameters[0],1.4823,1e-6);near(identified.parameters[1],.371,1e-5);
    // Independent nonlinear optimizer check, including an active bound.
    const auto rosen=[](const auto& x){return std::vector<double>{10*(x[1]-x[0]*x[0]),1-x[0]};};
    const auto optimum=numerics::boundedNonlinearLeastSquares(rosen,{-1.2,1},{-2,-1},{2,3});
    near(optimum.parameters[0],1,1e-5);near(optimum.parameters[1],1,1e-5);
    const auto bound=numerics::boundedNonlinearLeastSquares([](const auto& x){return std::vector<double>{x[0]-3};},{.5},{0},{1});
    near(bound.parameters[0],1);check(bound.converged);
    std::cout<<"Periodic interpolation, Bayes and causality tests passed\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
