#include "commands.hpp"
#include <dynsyslab/core/AcousticReceiver.hpp>
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/System.hpp>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace c=dynsyslab::core;
namespace fs=std::filesystem;
using Json=c::json::Value;
void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}
void close(double a,double b){require(std::abs(a-b)<1e-12+1e-12*std::abs(b),"Numeric fixture mismatch");}
void signalFixture(const fs::path& dir) {
    c::SystemConfiguration cfg;cfg.identifier="rpe";cfg.preset="air_water_20c";
    const auto system=c::createSystem(cfg);const double T=system->coordinates().secondsPerModelTime;
    Json::Object pars,groups;
    for(const auto& entry:system->configuration().parameters)pars[entry.first]=entry.second;
    for(const auto& entry:system->configuration().dimensionlessGroups)groups[entry.first]=entry.second;
    c::NpzRecord data;
    data.arrays["time"]=c::NpyArray::fromDoubles({5},{0.,1.,2.,3.,4.});
    const std::vector<double> state{1.2,0.,1.1,-1.,1.,1.,1.3,0.,1.15,-1.};
    data.arrays["values"]=c::NpyArray::fromDoubles({5,2},state);
    data.arrays["metadata_json"]=c::NpyArray::fromString(c::json::dump(Json::Object{
        {"system","rpe"},{"forcing_type","free"},{"integration_form","dimensionless"},{"system_preset",cfg.preset},
        {"time_scale_to_seconds",T},{"resolved_system_parameters",pars},{"dimensionless_groups",groups}}));
    c::writeNpzAtomicNew(dir/"state.npz",data);
    const double far=.005+500*T;  // Far channel lags by half a native time unit.
    const Json::Array positions{Json::Array{0.,far,0.},Json::Array{.005,0.,0.}};
    std::ofstream(dir/"geometry.json")<<c::json::dump(Json::Object{{"schema","dynsyslab.receiver-sidecar/v1"},
        {"stride",1.},{"sound_speed_m_s",1000.},{"positions_m",positions}});
    receiverSignals(dir/"state.npz",dir/"geometry.json",dir/"signals.npz");
    const auto signals=c::readNpz(dir/"signals.npz");
    auto windowed=c::json::load((dir/"geometry.json").wstring());
    windowed.object()["window_tau"]=Json::Array{1.5,3.};
    std::ofstream(dir/"window.json")<<c::json::dump(windowed);
    receiverSignals(dir/"state.npz",dir/"window.json",dir/"windowed.npz");
    const auto segment=c::readNpz(dir/"windowed.npz");
    require(segment.at("time").doubles()==std::vector<double>({2.,3.}),"Window changed the source clock");
    const auto clipped=segment.at("monopole_pressure").doubles(),full=signals.at("monopole_pressure").doubles();
    require(clipped==std::vector<double>(full.begin()+2,full.begin()+6),"Window lost earlier emission history");
    require(signals.at("time").doubles()==std::vector<double>({1.,2.,3.,4.}),"Missing propagation warmup");
    require(signals.at("ring_bounds").doubles()==std::vector<double>({0.,3.}),"Wrong radius extrema");
    const auto near=c::rpeReceiverPressures(*system,{{.005,0.,0.},1000.});
    const auto distant=c::rpeReceiverPressures(*system,{{0.,far,0.},1000.});
    for(const char* key:{"monopole_pressure","spherical_near_field_pressure","dipole_induction_pressure","dipole_radiation_pressure"}){
        const auto actual=signals.at(key).doubles(),past=signals.at(std::string("near_at_")+key).doubles();
        for (std::size_t i=0;i<4;++i) {
            auto value=[&](const c::ObservableCatalog& map,std::size_t j){return map.at(key).evaluateReal(&state[2*j])[0];};
            close(actual[2*i],.5*(value(distant,i)+value(distant,i+1)));
            close(past[2*i],.5*(value(near,i)+value(near,i+1)));
            close(actual[2*i+1],value(near,i+1));close(past[2*i+1],actual[2*i+1]);
        }
    }
}
void ringFixture(const fs::path& dir) {
    c::NpzRecord ref,fit;
    ref.arrays["time"]=c::NpyArray::fromDoubles({4},{.5,1.,1.5,2.});
    ref.arrays["ring_bounds"]=c::NpyArray::fromDoubles({3},{0.,1.,2.});
    ref.arrays["state"]=c::NpyArray::fromDoubles({4,1},{2.,1.,1.,4.});
    ref.arrays["error_center"]=c::NpyArray::fromDoubles({1},{1.});
    fit.arrays["time"]=ref.at("time");
    fit.arrays["prediction"]=c::NpyArray::fromDoubles({4,1},{1.5,1.,2.,4.});
    c::writeNpzAtomicNew(dir/"ref.npz",ref);c::writeNpzAtomicNew(dir/"fit.npz",fit);
    ringErrors(dir/"ref.npz",dir/"fit.npz",dir/"rings.npz");
    const auto rings=c::readNpz(dir/"rings.npz");
    require(rings.at("partial").doubles()==std::vector<double>({1.,0.}),"Partial ring mislabeled");
    require(rings.at("count").doubles()==std::vector<double>({1.,2.}),"Ring boundary double-counted");
    const auto score=rings.at("relative_rms_percent").doubles();
    close(score[0],50.);require(std::isnan(score[1]),"Zero reference deviation must be undefined");
}
int main(){
    const auto dir=fs::temp_directory_path()/("receiver-test-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    try{
        if(!fs::create_directory(dir))throw std::runtime_error("Fixture already exists");
        signalFixture(dir);ringFixture(dir);fs::remove_all(dir);return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
