#include <dynsyslab/core/AnalyticHistories.hpp>
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/System.hpp>
#include <cmath>
#include <stdexcept>
namespace dynsyslab::workflow {
core::json::Value evaluateAnalyticHistoryQuery(const core::json::Value&q,const core::ODESystem&system){
    using J=core::json::Value;using A=J::Array;using O=J::Object;
    const auto real=[](const J&v){if(!std::holds_alternative<double>(v.storage())||!std::isfinite(v.number()))throw std::invalid_argument("Finite history number required");return v.number();};
    const auto&states=q.at("initial_states");const auto&times=q.at("time_offsets");
    if(!std::holds_alternative<A>(states.storage())||!std::holds_alternative<A>(times.storage()))throw std::invalid_argument("Initial states and time offsets must be arrays");
    const auto dimension=system.descriptor().stateDimension;std::vector<double> initial,offsets;
    for(const auto&row:states.array()){if(!std::holds_alternative<A>(row.storage())||row.array().size()!=dimension)throw std::invalid_argument("Analytic history IC dimension mismatch");for(const auto&v:row.array())initial.push_back(real(v));}
    for(const auto&t:times.array())offsets.push_back(real(t));
    const auto catalog=system.observableCatalog();std::vector<core::Observable> observables;A names;
    for(const auto&name:q.at("observables").array()){const auto&o=catalog.at(name.string());observables.push_back(o);names.emplace_back(name.string());}
    const auto path=std::filesystem::u8path(q.at("output_archive").string());
    if(!path.empty()&&std::filesystem::exists(path))throw std::invalid_argument("History archive exists; use a new path");
    const double t0=q.at("initial_time").isNull()?0:real(q.at("initial_time"));
    if(path.empty()&&initial.size()/dimension*offsets.size()>100000)throw std::invalid_argument("Large histories require output_archive");
    const auto data=core::sampleAnalyticHistories(system.analyticSolutions(),initial,dimension,offsets,observables,t0,q.at("family").string("auto"));
    A families;for(const auto&f:data.solutionFamilies)families.emplace_back(f);
    O metadata{{"schema","dynsyslab.analytic-histories/v1"},{"system",system.configuration().identifier},
        {"equations_fingerprint",core::equationsFingerprint(system.configuration())},{"request",q},
        {"solution_families",families},{"layout","point/time/channel; complex128; physical/native time offsets relative to initial_time"},
        {"source","Registered analytic trajectory and observable evaluators; no numerical time integration or fitted dynamics"}};
    O result{{"status","complete"},{"points",double(data.points)},{"times",double(data.times)},{"channels",double(data.channels)},
        {"observables",names},{"metadata",metadata}};
    if(!path.empty()){core::NpzRecord archive;
        archive.arrays["values"]=core::NpyArray::fromComplex({data.points,data.times,data.channels},data.values);
        archive.arrays["initial_states"]=core::NpyArray::fromDoubles({data.points,dimension},initial);
        archive.arrays["time_offsets"]=core::NpyArray::fromDoubles({data.times},offsets);
        archive.arrays["metadata_json"]=core::NpyArray::fromString(core::json::dump(metadata));
        core::writeNpzAtomicNew(path,archive);result["archive"]=path.u8string();
    }else{A values;for(auto z:data.values)values.emplace_back(A{z.real(),z.imag()});result["values"]=values;}
    return result;
}
}
