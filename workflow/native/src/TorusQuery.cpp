#include "TorusQuery.hpp"
#include <dynsyslab/core/TorusAnalysis.hpp>
#include <dynsyslab/math/TorusFourier.hpp>
#include <cmath>
#include <stdexcept>
namespace dynsyslab::workflow {
namespace {
using J=core::json::Value;using A=J::Array;using O=J::Object;
double number(const J&v){if(!std::holds_alternative<double>(v.storage())||!std::isfinite(v.number()))throw std::invalid_argument("Finite number required");return v.number();}
int integer(const J&v,int minimum,int maximum){const double x=number(v);if(x<minimum||x>maximum||x!=std::floor(x))throw std::invalid_argument("Integer argument out of range");return int(x);}
template<class T>A values(const std::vector<T>&v){A out;for(auto x:v)out.emplace_back(double(x));return out;}
}
core::json::Value evaluateTorusAnalysis(const J&q,const core::ODESystem&s){
    O out{{"status","not_applicable"},{"integrability_certified",false},{"kam_certified",false},
        {"qualification_scope","A registered analytic invariant-set evaluator plus numerical checks; no automatic proof of global integrability"}};
    core::AnalyticState initial;for(const auto&v:q.at("initial_state").array())initial.push_back(number(v));
    if(initial.size()!=s.descriptor().stateDimension)throw std::invalid_argument("Initial state dimension mismatch");
    const double t0=q.at("initial_time").isNull()?0:number(q.at("initial_time"));core::AnalyticTrajectory a;std::string reason;
    const auto selection=q.at("family").string("auto");
    for(const auto&f:s.analyticSolutions()){
        if(selection!="auto"&&selection!=f.descriptor.identifier)continue;
        auto candidate=f.instantiate(initial,t0);
        if(candidate&&(candidate.torus||(std::isfinite(candidate.period)&&candidate.period>0))){a=std::move(candidate);out["family"]=f.descriptor.identifier;break;}
        reason+=f.descriptor.identifier+": "+(candidate?"No regular compact torus/period supplied":candidate.reason)+"; ";
    }
    if(!a){out["reason"]=reason.empty()?"No applicable analytic family is registered":reason;return out;}
    const auto dimensions=a.angularFrequencies.empty()?1:a.angularFrequencies.size();std::vector<std::size_t>shape;
    if(!q.at("shape").isNull()){if(!std::holds_alternative<A>(q.at("shape").storage()))throw std::invalid_argument("shape must be an array");
        for(const auto&v:q.at("shape").array())shape.push_back(std::size_t(integer(v,3,1000000)));}
    else shape.assign(dimensions,dimensions==1?256:32);
    std::vector<core::Observable>obs;A names;const auto catalog=s.observableCatalog();
    if(q.at("observables").isNull()){for(const auto&n:s.descriptor().stateNames){obs.push_back(catalog.at(n));names.emplace_back(n);}}
    else for(const auto&name:q.at("observables").array()){obs.push_back(catalog.at(name.string()));names.push_back(name);}
    auto sampled=core::sampleInvariantTorus(s,a,obs,shape,t0);
    out["status"]="complete";out["torus_dimension"]=double(dimensions);out["angular_frequencies"]=values(sampled.angularFrequencies);
    out["shape"]=values(shape);out["sample_convention"]="uniform-product-angles-last-axis-fastest-no-endpoint";
    out["state_dimension"]=double(sampled.stateDimension);out["values"]=values(sampled.states);out["observables"]=names;
    out["maximum_invariance_residual"]=sampled.maximumInvarianceResidual;out["maximum_wrap_residual"]=sampled.maximumWrapResidual;
    out["sampled_checks_pass"]=sampled.maximumInvarianceResidual<1e-6&&sampled.maximumWrapResidual<1e-8;
    A measurement;for(auto v:sampled.observables)measurement.emplace_back(A{v.real(),v.imag()});out["measurements"]=measurement;
    std::vector<int>maximum;
    if(!q.at("maximum_harmonics").isNull())for(const auto&v:q.at("maximum_harmonics").array())maximum.push_back(integer(v,0,1000000));
    else for(auto n:shape)maximum.push_back(int((n-1)/2));
    A series;for(std::size_t j=0;j<obs.size();++j){std::vector<std::complex<double>>channel;
        for(std::size_t i=j;i<sampled.observables.size();i+=obs.size())channel.push_back(sampled.observables[i]);
        const auto f=math::torusFourierSeries(channel,shape,sampled.angularFrequencies,maximum);A coefficients;
        for(const auto&[k,c]:f.coefficients)coefficients.emplace_back(O{{"index",values(k)},{"coefficient",A{c.real(),c.imag()}},{"angular_frequency",math::harmonicFrequency(k,f.angularFrequencies)}});
        series.emplace_back(O{{"angular_frequencies",values(f.angularFrequencies)},{"coefficients",coefficients},{"mean_square",math::meanSquare(f)}});
    }out["torus_fourier_series"]=series;return out;
}
}
