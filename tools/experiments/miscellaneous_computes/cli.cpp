#include "hankel.hpp"
#include "rollout_cli.hpp"
#include "crossing_cli.hpp"
#include "local_rhs.hpp"
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/SindyDynamics.hpp>
#include <dynsyslab/core/CoordinateTransform.hpp>
#include <cmath>
#include <climits>
#include <filesystem>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

namespace core = dynsyslab::core;
namespace fs = std::filesystem;
namespace {
std::size_t natural(const std::string& value) {
    std::size_t used = 0;
    auto result = std::stoull(value,&used);
    if (value.empty() || value.front()=='-' || used!=value.size() || !result)
        throw std::invalid_argument("Expected a positive integer");
    return result;
}
void write(core::NpzRecord& record,const fs::path& path,const std::string& kind,const fs::path& input) {
    record.arrays["metadata_json"]=core::NpyArray::fromString(core::json::dump(core::json::Value::Object{
        {"schema","dynsyslab.miscellaneous-computes/"+kind+"/v1"},{"input",fs::absolute(input).generic_u8string()}}));
    core::writeNpzAtomicNew(path,record);
}
void gram(const fs::path& path,std::size_t taps,const fs::path& dest) {
    auto input=core::readNpz(path);
    const auto& arr=input.at("values");
    if(arr.shape.size()!=2 || !arr.shape[1] || taps>arr.shape[0])
        throw std::invalid_argument("values must be samples by channels and include a complete window");
    auto d=taps*arr.shape[1];
    if(d/taps!=arr.shape[1] || d>std::vector<double>().max_size()/d)
        throw std::length_error("Gram matrix exceeds addressable memory");
    std::vector<double> G(d*d);
    if(mc_hankel_gram_v1(arr.doubles().data(),arr.shape[0],arr.shape[1],taps,G.data()))
        throw std::invalid_argument("Invalid data or insufficient workspace for Hankel covariance");
    core::NpzRecord result; result.arrays["gram"]=core::NpyArray::fromDoubles({d,d},G);
    write(result,dest,"gram",path);
}
void integral(const fs::path& path,double duration,const fs::path& dest) {
    const auto input=core::readNpz(path);
    auto t=input.at("time").doubles(),y=input.at("values").doubles();
    const auto& shape=input.at("values").shape;
    if(input.at("time").shape!=std::vector<std::size_t>{t.size()} ||
       (shape!=std::vector<std::size_t>{t.size()} && shape!=std::vector<std::size_t>{t.size(),1}))
        throw std::invalid_argument("Window integral requires aligned scalar observations");
    auto q=core::windowIntegralCoordinates(core::TimeSeries(std::move(t),std::move(y),1),duration);
    const auto n=q.sampleCount();core::NpzRecord result;
    result.arrays["time"]=core::NpyArray::fromDoubles({n},{q.timeData(),q.timeData()+n});
    result.arrays["values"]=core::NpyArray::fromDoubles({n,2},{q.valuesData(),q.valuesData()+2*n});
    result.arrays["duration"]=core::NpyArray::fromDoubles({1},{duration});
    write(result,dest,"window-integral",path);
}
void fit(const fs::path& path,std::size_t degree,const std::string& threshold,const fs::path& dest) {
    auto input=core::readNpz(path);
    const auto& b=input.at("b");
    if(b.shape.size()!=2 || !b.shape[1] || b.shape!=input.at("db").shape)
        throw std::invalid_argument("Expected matching two-dimensional b and db arrays");
    auto t=input.at("time").doubles();
    if(t.size()!=b.shape[0])throw std::invalid_argument("Clock length mismatch");
    core::TimeSeries coordinates(t,b.doubles(),b.shape[1]),rates(t,input.at("db").doubles(),b.shape[1]);
    core::SindyDynamicsOptions opts; std::size_t used=0;opts.threshold=std::stod(threshold,&used);
    if(used!=threshold.size() || !std::isfinite(opts.threshold) || opts.threshold<0)
        throw std::invalid_argument("Expected a finite nonnegative threshold");
    auto model=core::fitSindyDynamics(coordinates,rates,core::polynomialDynamicsLibrary(b.shape[1],degree),opts);
    core::NpzRecord result; result.arrays["model"]=core::NpyArray::fromString(core::json::dump(model.definition()));
    write(result,dest,"sindy",path);
}
void roll(const fs::path& modelPath,const fs::path& queryPath,const fs::path& dest) {
    auto saved=core::readNpz(modelPath),query=core::readNpz(queryPath);
    auto model=core::SindyDynamicsModel::fromDefinition(core::json::parse(saved.at("model").string()));
    auto t=query.at("time").doubles(),b=query.at("initial").doubles();
    if(t.size()<2 || b.size()!=model.dimension())throw std::invalid_argument("Invalid rollout dimensions");
    core::TimeSeries clock(t,std::vector<double>(t.size()),1);
    if(!clock.uniformTimeSteps())throw std::invalid_argument("The native solver gate expects a uniform output clock");
    core::IntegrationOptions opts;opts.method=core::IntegrationMethod::dormandPrince45;
    opts.initialTime=t.front();opts.finalTime=t.back();opts.outputStep=clock.deltaTime();
    opts.relativeTolerance=1.e-9;opts.absoluteTolerance=1.e-11;
    auto solution=core::integrateSindyDynamics(model,b,opts);
    core::NpzRecord result;
    result.arrays["time"]=core::NpyArray::fromDoubles({solution.time.size()},solution.time);
    result.arrays["b"]=core::NpyArray::fromDoubles({solution.time.size(),b.size()},solution.values);
    result.arrays["model_source"]=core::NpyArray::fromString(fs::absolute(modelPath).generic_u8string());
    write(result,dest,"rollout",queryPath);
}
void localRates(const fs::path& fieldPath,const fs::path& queryPath,const fs::path& dest) {
    const auto field=core::readNpz(fieldPath),query=core::readNpz(queryPath);
    if(field.at("schema").string()!="dynsyslab.miscellaneous-computes/local-field/v1")
        throw std::invalid_argument("Unsupported local field schema");
    const auto shape=field.at("center").shape,states=query.at("b").shape;
    if(shape.size()!=2 || states.size()!=2 || states[1]!=shape[1] || !states[0] ||
       states[0]>INT_MAX || shape[0]>INT_MAX || shape[1]>64)
        throw std::invalid_argument("Invalid local field or query dimensions");
    const auto rows=shape[0],n=shape[1];
    auto checked=[&](const char* key,std::vector<std::size_t> expected) {
        const auto& array=field.at(key);
        if(array.shape!=expected)throw std::invalid_argument(std::string("Invalid array shape for ")+key);
        return array.doubles();
    };
    auto centers=checked("center",shape),means=checked("mean",shape),p=checked("P",{rows,n,2}),
         q=checked("Q",{rows,n,2}),coef=checked("coef",{rows,6,n}),scale=checked("scale",{n});
    const double k=checked("neighbors",{}).front(),d=checked("decay",{}).front();
    if(!(k>=1 && k<=64 && k==std::floor(k)))throw std::invalid_argument("Invalid neighborhood size");
    auto x=query.at("b").doubles();std::vector<double> derivatives(x.size());
    const auto finite=[](const auto& a){for(double v:a)if(!std::isfinite(v))return false;return true;};
    if(!finite(x)||!finite(centers)||!finite(means)||!finite(p)||!finite(q)||!finite(coef)||!finite(scale))
        throw std::invalid_argument("Nonfinite local model or query");
    for(double v:scale)if(v<=0.)throw std::invalid_argument("Invalid local normalization");
    std::unique_ptr<void,decltype(&mc_local_drop_v1)> engine(mc_local_field_v1(static_cast<int>(n),
        static_cast<int>(rows),static_cast<int>(k),d,centers.data(),means.data(),p.data(),q.data(),
        coef.data(),scale.data()),mc_local_drop_v1);
    if(!engine)throw std::invalid_argument("Unable to construct local field");
    mc_local_batch_v1(engine.get(),static_cast<int>(states[0]),x.data(),derivatives.data());
    core::NpzRecord result;
    result.arrays["db"]=core::NpyArray::fromDoubles(states,derivatives);
    result.arrays["model_source"]=core::NpyArray::fromString(fs::absolute(fieldPath).generic_u8string());
    write(result,dest,"local-evaluation",queryPath);
}
}
int main(int argc,char** argv) {
    try {
        const std::string command=argc>1?argv[1]:"";
        if(command=="gram" && argc==5)gram(argv[2],natural(argv[3]),argv[4]);
        else if(command=="crossings" && argc==4)crossingTable(argv[2],argv[3]);
        else if(command=="integral" && argc==5){
            const std::string value=argv[3];std::size_t used=0;double duration=std::stod(value,&used);
            if(used!=value.size() || !std::isfinite(duration) || duration<=0.)throw std::invalid_argument("Invalid window duration");
            integral(argv[2],duration,argv[4]);
        }
        else if(command=="fit" && argc==6)fit(argv[2],natural(argv[3]),argv[4],argv[5]);
        else if(command=="rollout" && argc==5)roll(argv[2],argv[3],argv[4]);
        else if(command=="rk4" && argc==6){
            std::size_t used=0;const std::string value=argv[4];
            const double step=std::stod(value,&used);
            if(used!=value.size())throw std::invalid_argument("Invalid maximum step");
            fixedRollout(argv[2],argv[3],step,argv[5]);
        }
        else if(command=="local-eval" && argc==5)localRates(argv[2],argv[3],argv[4]);
        else {
            std::cout<<"misc_computes gram signal.npz taps out.npz\n"
                     <<"misc_computes crossings plane.npz out.npz\n"
                     <<"misc_computes integral signal.npz duration out.npz\n"
                     <<"misc_computes rk4 model.npz initial-and-clock.npz max-step out.npz\n"
                     <<"misc_computes fit coordinates.npz degree threshold model.npz\n"
                     <<"misc_computes rollout model.npz initial-and-clock.npz out.npz\n"
                     <<"misc_computes local-eval local-field.npz query.npz rates.npz\n";
            return command=="--help"?0:2;
        }
    } catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
    return 0;
}
