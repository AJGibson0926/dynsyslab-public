// Measurement-only dense paired input fitting. Echo prehistory is referenced,
// not copied into thousands of overlapping archives. No state access here.
#include "echo_pulse.hpp"
#include <cmath>
#include <stdexcept>
namespace echo_pulse {
std::map<std::string,ModeSamples> denseSamples(const fs::path& root,const fs::path& out,
    const std::string& window,const Frozen& f,const core::KoopmanModel& bank) {
    std::map<std::string,ModeSamples> result;Json::Array diagnostics;
    const auto folder=root/"archives/experiments/echo-input-regression-20261006-v1";
    const auto manifest=readJson(out/"dense-acquisition.json");
    for(const auto& meta:manifest.at("records").array()) {
        const auto id=meta.at("source").string(),key=meta.at("split").string()+"-"+meta.at("regime").string();
        const auto dest=out/(window+"-dense-"+id+".samples.npz");
        if(result[key].empty())result[key].resize(3);
        core::NpzRecord samples;
        if(fs::exists(dest))samples=core::readNpz(dest);
        else {
            const auto a=core::readNpz(out/("dense-"+id+".pulse.npz"));
            if(a.contains("states"))throw std::runtime_error("dense analysis rejects state-bearing input");
            const auto source=readMeasurements(folder/(id+".measurements.npz"));
            const auto indices=a.at("source_indices").doubles(),tails=a.at("echo_tails").doubles();
            std::vector<Complex> values,sensitivity;std::vector<double> coordinates;
            core::CausalHistoryProjector projector(f.basis);std::size_t seen=0;
            for(std::size_t k=0;k<indices.size();++k) {
                const auto index=std::size_t(indices[k]);while(seen<=index)projector.append(source.echo.at(seen++));
                const auto xi=projector.projectLatest();const auto z=bank.evaluateCoordinates(xi->data());
                std::vector<Complex> after[2];
                for(std::size_t b=0;b<2;++b) {
                    auto continuation=projector;
                    for(std::size_t s=1;s<=f.stride;++s)continuation.append(tails[(k*2+b)*(f.stride+1)+s]);
                    const auto next=continuation.projectLatest();after[b]=bank.evaluateCoordinates(next->data());
                }
                coordinates.push_back((*xi)[0]);coordinates.push_back((*xi)[1]);
                for(std::size_t j=0;j<3;++j){values.push_back(z[j]);sensitivity.push_back((after[0][j]-after[1][j])/(2*meta.at("delta").number()));}
            }
            samples.arrays["psi"]=core::NpyArray::fromComplex({indices.size(),3},values);
            samples.arrays["sensitivity"]=core::NpyArray::fromComplex({indices.size(),3},sensitivity);
            samples.arrays["xi"]=core::NpyArray::fromDoubles({indices.size(),2},coordinates);
            samples.arrays["metadata_json"]=core::NpyArray::fromString(core::json::dump(meta));
            core::writeNpzAtomicNew(dest,samples);
        }
        const auto z=samples.at("psi").complexDoubles(),b=samples.at("sensitivity").complexDoubles();
        const auto x=samples.at("xi").doubles();double cancellation=0;
        for(std::size_t k=0;k<z.size()/3;++k)for(std::size_t j=0;j<3;++j) {
            const auto mu=std::exp(bank.continuousEigenvalues()[j]*.01),next=mu*z[3*k+j]+b[3*k+j];
            cancellation=std::max(cancellation,std::abs((next-mu*z[3*k+j])-b[3*k+j])/std::max(1e-300,std::abs(b[3*k+j])));
            result[key][j].push_back({{x[2*k],x[2*k+1]},z[3*k+j],next,1.});
        }
        if(cancellation>1e-5)throw std::runtime_error("dense sensitivity surrogate precision loss");
        diagnostics.emplace_back(Json::Object{{"source",id},{"points",double(z.size()/3)},{"surrogate_relative_cancellation",cancellation}});
    }
    const auto check=out/(window+"-dense-projection.json");if(!fs::exists(check))saveJson(check,Json::Object{{"sources",diagnostics}});
    return result;
}
}
