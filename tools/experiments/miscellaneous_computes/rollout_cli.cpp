#include "rollout_cli.hpp"
#include "rollout_kernel.hpp"
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/SindyDynamics.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace c = dynsyslab::core;
namespace fs = std::filesystem;
namespace {
int evaluate(void* ctx,const double* x,double* out) {
    const auto rates=static_cast<c::SindyDynamicsModel*>(ctx)->evaluate(x);
    std::copy(rates.begin(),rates.end(),out);
    return 0;
}
}
void fixedRollout(const fs::path& modelPath,const fs::path& queryPath,double step,const fs::path& dest) {
    const auto stored=c::readNpz(modelPath),query=c::readNpz(queryPath);
    auto field=c::SindyDynamicsModel::fromDefinition(c::json::parse(stored.at("model").string()));
    auto x=query.at("initial").doubles(),time=query.at("time").doubles();
    const auto n=field.dimension();
    if (!n || query.at("initial").shape!=std::vector<std::size_t>{n} ||
        query.at("time").shape!=std::vector<std::size_t>{time.size()} || time.empty() ||
        time.size()>std::vector<double>().max_size()/n)
        throw std::invalid_argument("Invalid initial state or output clock");
    std::vector<double> limits;
    if (query.contains("limits")) {
        limits=query.at("limits").doubles();
        if (limits.size()!=n) throw std::invalid_argument("Expected one magnitude limit per state");
    }
    std::size_t budget=100000000;
    if (query.contains("max_steps")) {
        auto values=query.at("max_steps").doubles();
        if (values.size()!=1 || !(values[0]>=1 && values[0]<=9007199254740991. &&
              values[0]<static_cast<double>(std::numeric_limits<std::size_t>::max())) || values[0]!=std::floor(values[0]))
            throw std::invalid_argument("Invalid max_steps");
        budget=static_cast<std::size_t>(values[0]);
    }
    std::vector<double> rows(n*time.size()); std::size_t written=0,evaluations=0;
    int code=mc_rk4_v1(evaluate,&field,n,x.data(),time.data(),time.size(),step,budget,
                        limits.empty()?nullptr:limits.data(),rows.data(),&written,&evaluations);
    if (code>=4) throw std::runtime_error("Invalid or unallocatable fixed-step request");
    rows.resize(written*n);time.resize(written);
    c::NpzRecord result;
    result.arrays["time"]=c::NpyArray::fromDoubles({written},time);
    result.arrays["b"]=c::NpyArray::fromDoubles({written,n},rows);
    result.arrays["metadata_json"]=c::NpyArray::fromString(c::json::dump(c::json::Value::Object{
        {"schema","dynsyslab.miscellaneous-computes/fixed-rollout/v1"},
        {"model",fs::absolute(modelPath).generic_u8string()},{"query",fs::absolute(queryPath).generic_u8string()},
        {"status",static_cast<double>(code)},{"complete",code==0},{"evaluations",static_cast<double>(evaluations)},{"max_step",step}}));
    c::writeNpzAtomicNew(dest,result);
}
