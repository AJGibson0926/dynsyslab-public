#include <dynsyslab/core/AnalyticHistories.hpp>
#include <cmath>
#include <stdexcept>
// Standard delay-observable construction g o Phi^t. Analytic evaluators and
// observable definitions own their equation-specific references. See Arbabi
// and Mezic, https://doi.org/10.1137/16M1105230 for delay Koopman dictionaries.
namespace dynsyslab::core {
AnalyticHistories sampleAnalyticHistories(const std::vector<AnalyticSolutionFamily>&families,
    const std::vector<double>&initial,std::size_t dimension,const std::vector<double>&offsets,
    const std::vector<Observable>&observables,double t0,const std::string&selection){
    if(!dimension||initial.empty()||initial.size()%dimension||offsets.empty()||observables.empty()
        ||!std::isfinite(t0))throw std::invalid_argument("Analytic history dimensions/clock are invalid");
    AnalyticHistories out;out.points=initial.size()/dimension;out.times=offsets.size();
    if(out.points>1000000||out.times>1000000)throw std::invalid_argument("Analytic history point/time limit exceeded");
    for(const auto&o:observables){if(o.stateDimension()!=dimension)throw std::invalid_argument("History observable state dimension mismatch");out.channels+=o.outputDimension();}
    if(out.channels>1024||out.times>50000000/out.channels||out.points>50000000/(out.times*out.channels))
        throw std::invalid_argument("Analytic history block exceeds fifty million complex entries; supply smaller IC blocks");
    for(double v:initial)if(!std::isfinite(v))throw std::invalid_argument("Finite initial states required");
    for(double t:offsets)if(!std::isfinite(t)||!std::isfinite(t0+t))throw std::invalid_argument("Finite history times required");
    out.values.reserve(out.points*out.times*out.channels);out.solutionFamilies.reserve(out.points);
    for(std::size_t i=0;i<out.points;++i){
        const AnalyticState state(initial.begin()+i*dimension,initial.begin()+(i+1)*dimension);
        AnalyticTrajectory orbit;std::string reason;
        for(const auto&family:families){if(selection!="auto"&&family.descriptor.identifier!=selection)continue;
            orbit=family.instantiate(state,t0);if(orbit){out.solutionFamilies.push_back(family.descriptor.identifier);break;}
            reason+=orbit.reason+"; ";}
        if(!orbit)throw std::invalid_argument("No applicable analytic history at IC "+std::to_string(i)+": "+reason);
        for(double t:offsets){const auto x=orbit.evaluate(t0+t);
            if(x.size()!=dimension)throw std::runtime_error("Analytic evaluator returned wrong dimension");
            for(double v:x)if(!std::isfinite(v))throw std::domain_error("Analytic history left the evaluator domain");
            for(const auto&o:observables){const auto values=o.evaluate(x.data());
                for(auto v:values){if(!std::isfinite(v.real())||!std::isfinite(v.imag()))throw std::domain_error("Nonfinite analytic history observable");out.values.push_back(v);}}
        }
    }return out;
}
}
