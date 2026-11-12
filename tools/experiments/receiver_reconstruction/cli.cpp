#include "lookup.hpp"
#include "commands.hpp"
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Json.hpp>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace c=dynsyslab::core;
namespace fs=std::filesystem;
namespace {
void estimate(const fs::path& trainPath,const fs::path& queryPath,std::size_t k,const fs::path& destination) {
    const auto train=c::readNpz(trainPath), query=c::readNpz(queryPath);
    const auto& shape=train.at("features").shape;
    const auto& states=train.at("state").shape;
    const auto& samples=query.at("features").shape;
    if(shape.size()!=2 || states.size()!=2 || samples.size()!=2 || shape[0]!=states[0] ||
       samples[1]!=shape[1] || !samples[0] || !k || k>shape[0] ||
       train.at("scales").shape!=std::vector<std::size_t>{shape[1]} ||
       query.at("time").shape!=std::vector<std::size_t>{samples[0]})
        throw std::runtime_error("Invalid training or query matrix dimensions");
    receivers::Lookup index(train.at("features").doubles(),train.at("state").doubles(),
                            train.at("scales").doubles(),states[1]);
    const auto x=query.at("features").doubles();
    std::vector<double> predictions,weights,nearest;std::vector<std::int64_t> rows;
    for(std::size_t i=0;i<samples[0];++i) {
        const auto est=index.predict(x.data()+i*samples[1],k);
        predictions.insert(predictions.end(),est.values.begin(),est.values.end());
        weights.insert(weights.end(),est.weights.begin(),est.weights.end());
        for(const auto& neighbor:est.neighbors)rows.push_back(static_cast<std::int64_t>(neighbor.row));
        nearest.push_back(std::sqrt(est.neighbors.front().squaredDistance));
    }
    c::NpzRecord out;
    out.arrays["time"]=query.at("time");out.arrays["scales"]=train.at("scales");
    out.arrays["prediction"]=c::NpyArray::fromDoubles({samples[0],states[1]},predictions);
    out.arrays["neighbor_rows"]=c::NpyArray::fromInt64s({samples[0],k},rows);
    out.arrays["weights"]=c::NpyArray::fromDoubles({samples[0],k},weights);
    out.arrays["nearest_distance"]=c::NpyArray::fromDoubles({samples[0]},nearest);
    out.arrays["metadata_json"]=c::NpyArray::fromString(c::json::dump(c::json::Value::Object{
        {"schema","dynsyslab.neighbor-estimate/v1"},{"k",double(k)},{"train",fs::absolute(trainPath).generic_u8string()},
        {"query",fs::absolute(queryPath).generic_u8string()},{"metric","Euclidean on asinh(feature/scale)"},
        {"weights","inverse square, distance floor 1e-12"}}));
    c::writeNpzAtomicNew(destination,out);
}
}
int main(int argc,char** argv) {
    try {
        if(argc==5 && std::string(argv[1])=="signals")receiverSignals(argv[2],argv[3],argv[4]);
        else if(argc==5 && std::string(argv[1])=="rings")ringErrors(argv[2],argv[3],argv[4]);
        else if(argc==6 && std::string(argv[1])=="estimate") {
            const std::string count=argv[4];std::size_t used=0;
            const auto k=std::stoull(count,&used);
            if(used!=count.size() || count[0]=='-')throw std::runtime_error("Invalid neighbor count");
            estimate(argv[2],argv[3],k,argv[5]);
        } else throw std::runtime_error("Usage: receiver_sidecar signals RPE.npz geometry.json out.npz | estimate train.npz query.npz K out.npz | rings reference.npz estimates.npz out.npz");
        return 0;
    } catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
