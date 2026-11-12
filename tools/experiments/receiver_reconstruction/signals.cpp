#include <dynsyslab/core/AcousticReceiver.hpp>
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/System.hpp>
#include <dynsyslab/core/TimeSeries.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <limits>
#include <stdexcept>

namespace c=dynsyslab::core;
namespace fs=std::filesystem;
using Json=c::json::Value;
namespace {
const std::array<const char*,4> names{"monopole_pressure","spherical_near_field_pressure",
    "dipole_induction_pressure","dipole_radiation_pressure"};
struct Channel {
    double distance, delay;
    c::ObservableCatalog terms;
    std::size_t left=0;
};
std::vector<double> maxima(const std::vector<double>& t,const std::vector<double>& x) {
    std::vector<double> result;
    if (x[1]==0. && x[3]<0.) result.push_back(t.front());
    for (std::size_t i=1;i<t.size();++i) {
        const double a=x[2*i-1], b=x[2*i+1];
        if (a>0 && b<=0) result.push_back(t[i-1]+a/(a-b)*(t[i]-t[i-1]));
    }
    return result;
}
double interp(const c::Observable& term, const double* a, const double* b, double f) {
    const double y=term.evaluateReal(a).front();
    return f==0 ? y : y+f*(term.evaluateReal(b).front()-y);
}
}

void receiverSignals(const fs::path& path,const fs::path& config,const fs::path& output) {
    const auto job=c::json::load(config.wstring());
    if(job.at("schema").string()!="dynsyslab.receiver-sidecar/v1")throw std::runtime_error("Unknown receiver configuration");
    const double step=job.at("stride").number(),speed=job.at("sound_speed_m_s").number();
    if(!(step>=1 && step<=1e9) || std::floor(step)!=step)throw std::runtime_error("Expected positive integer stride");
    const auto stride=static_cast<std::size_t>(step);
    double begin=-std::numeric_limits<double>::infinity(),end=-begin;
    const auto& window=job.at("window_tau");
    if(!window.isNull()) {
        if(window.array().size()!=2)throw std::runtime_error("window_tau requires two endpoints");
        begin=window.array()[0].number(std::numeric_limits<double>::quiet_NaN());
        end=window.array()[1].number(std::numeric_limits<double>::quiet_NaN());
        if(!std::isfinite(begin) || !std::isfinite(end) || end<begin)
            throw std::runtime_error("Invalid emission time window");
    }
    const auto record=c::readNpz(path);
    const auto meta=c::json::parse(record.metadataJson());
    if(meta.at("system").string()!="rpe" || meta.at("integration_form").string()!="dimensionless" ||
       meta.at("archive_coordinate_system").string("dimensionless")!="dimensionless" ||
       meta.at("forcing_type").string()!="free" || !meta.at("coordinate_provider").isNull())
        throw std::runtime_error("Expected an unforced physical RPE trajectory in native dimensionless coordinates");
    const auto t=record.at("time").doubles(),x=record.at("values").doubles();
    if(t.size()<2 || record.at("time").shape.size()!=1 ||
       record.at("values").shape!=std::vector<std::size_t>{t.size(),2})
        throw std::runtime_error("Invalid physical state archive shape");
    (void)c::TimeSeries::view(t,x.data(),t.size(),2);
    c::SystemConfiguration sc;sc.identifier="rpe";sc.preset=meta.at("system_preset").string();
    const auto& pars=meta.at("resolved_system_parameters");
    if(pars.object().empty())throw std::runtime_error("Missing resolved physical parameters");
    for(const auto& item:pars.object())sc.parameters[item.first]=item.second.number(std::numeric_limits<double>::quiet_NaN());
    const auto system=c::createSystem(sc);
    const auto& cfg=system->configuration();
    const double T=system->coordinates().secondsPerModelTime,R0=cfg.parameters.at("equilibrium_radius_m");
    if(std::abs(meta.at("time_scale_to_seconds").number()-T)>1e-9*T)
        throw std::runtime_error("Archive clock does not match the reconstructed RPE");
    std::vector<Channel> channels;
    for(const auto& value:job.at("positions_m").array()) {
        const auto& a=value.array();
        if(a.size()!=3)throw std::runtime_error("Expected three position components");
        const double bad=std::numeric_limits<double>::quiet_NaN();
        const c::AcousticReceiver receiver{{a[0].number(bad),a[1].number(bad),a[2].number(bad)},speed};
        channels.push_back({receiver.rangeMetres(),receiver.travelTimeSeconds()/T,c::rpeReceiverPressures(*system,receiver)});
    }
    if(channels.empty())throw std::runtime_error("No receivers configured");
    const auto nearest=std::min_element(channels.begin(),channels.end(),[](const Channel& a,const Channel& b){return a.distance<b.distance;});
    const double nearDelay=nearest->delay;
    for(std::size_t i=0;i<t.size();++i)
        if(!(x[2*i]>0 && x[2*i]*R0<nearest->distance) || !std::isfinite(x[2*i+1]))
            throw std::runtime_error("Invalid physical sample or receiver inside bubble");
    std::vector<double> ranges,lag;
    for(const auto& chan:channels){ranges.push_back(chan.distance);lag.push_back(chan.delay-nearDelay);}
    const double maxLag=*std::max_element(lag.begin(),lag.end());
    std::vector<double> times,states;
    std::array<std::vector<double>,4> sensor,nearAtEmission;
    for(std::size_t i=0;i<t.size();i+=stride) {
        if(t[i]>end)break;
        if(t[i]<begin || t[i]-maxLag<t.front())continue;
        times.push_back(t[i]);states.insert(states.end(),{x[2*i],x[2*i+1]});
        for(std::size_t j=0;j<channels.size();++j) {
            auto& ch=channels[j];const double ts=t[i]-lag[j];
            while(ch.left+1<t.size() && t[ch.left+1]<=ts)++ch.left;
            const auto right=std::min(ch.left+1,t.size()-1);
            const double fraction=right==ch.left?0.:(ts-t[ch.left])/(t[right]-t[ch.left]);
            const auto* a=x.data()+2*ch.left;const auto* b=x.data()+2*right;
            for(std::size_t k=0;k<names.size();++k) {
                const double pressure=interp(ch.terms.at(names[k]),a,b,fraction);
                sensor[k].push_back(pressure);
                nearAtEmission[k].push_back(j==static_cast<std::size_t>(nearest-channels.begin())?
                    pressure:interp(nearest->terms.at(names[k]),a,b,fraction));
            }
        }
    }
    if(times.empty())throw std::runtime_error("No samples in the requested window after propagation warmup");
    c::NpzRecord out;
    auto put=[&](const std::string& key,const std::vector<double>& v,std::vector<std::size_t> shape){
        out.arrays[key]=c::NpyArray::fromDoubles(shape,v);
    };
    const auto n=times.size(),d=channels.size();
    const auto bounds=maxima(t,x);
    put("time",times,{n});put("state",states,{n,2});put("ring_bounds",bounds,{bounds.size()});
    put("error_center",{1.,0.},{2});
    put("ranges_m",ranges,{d});put("lag_tau",lag,{d});
    for(std::size_t k=0;k<names.size();++k) {
        put(names[k],sensor[k],{n,d});put(std::string("near_at_")+names[k],nearAtEmission[k],{n,d});
    }
    out.arrays["metadata_json"]=c::NpyArray::fromString(c::json::dump(Json::Object{
        {"schema","dynsyslab.receiver-signals/v1"},{"source",fs::absolute(path).generic_u8string()},
        {"clock","near receiver emission, native tau"},{"units","Pa"},{"time_scale_s",T},
        {"dimensionless_groups",meta.at("dimensionless_groups")},{"parameters",pars},{"receiver_config",job}}));
    c::writeNpzAtomicNew(output,out);
}
