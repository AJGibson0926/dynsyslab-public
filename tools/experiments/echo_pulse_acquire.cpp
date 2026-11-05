// ONLY acquisition functions may access archived plant states or a simulator.
// The returned archive contains echo histories and echo continuations only.
#include "echo_pulse.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
namespace echo_pulse {
Json readJson(const fs::path& p){std::ifstream f(p);std::stringstream s;s<<f.rdbuf();return core::json::parse(s.str());}
std::vector<Preparation> preparations() {
    std::vector<Preparation> out;
    auto add=[&](const std::string& source,const std::string& regime,const std::string& split,
                 const std::vector<double>& times) {
        for(std::size_t j=0;j<times.size();++j)out.push_back({source+"-p"+std::to_string(j),source,regime,split,
            std::round(times[j]*100)/100,regime=="nonlinear"?.1:.01,j==0,j==0&&source.back()=='1'});
    };
    for(int n:{1,3,5})add("free-train-"+std::to_string(n),"free","train",{8.13,16.37,35.69,50.11});
    for(int n:{1,3,5})add("free-validation-"+std::to_string(n),"free","heldout",{12.23,42.47});
    for(const std::string regime:{"mild","nonlinear"}) {
        for(int n:{1,3,5})add("local-"+regime+"-"+std::to_string(n),regime,"train",{20.13,24.37,28.69,33.11});
        for(int n:{2,4,6,7,8})add("local-"+regime+"-"+std::to_string(n),regime,"heldout",{21.23,30.47});
        std::vector<double> times;for(int j=0;j<12;++j)times.push_back(30+j/(12*(regime=="mild"?.93:.83)));
        add("target-"+regime,regime,"target-diagnostic",times);
    }
    return out;
}
Json acquire(const fs::path& root,const fs::path& output,const Preparation& p,const Frozen& f,
             const core::SystemConfiguration& config) {
    const auto dest=output/(p.id+".pulse.npz");
    if(fs::exists(dest))return core::json::parse(core::readNpz(dest).at("metadata_json").string());
    const auto folder=root/"archives/experiments/echo-input-regression-20261006-v1";
    const auto rec=readMeasurements(folder/(p.source+".measurements.npz"));
    const auto truth=core::readNpz(folder/(p.source+".truth.npz"));
    if(rec.time!=truth.at("time").doubles())throw std::runtime_error("acquisition source timing mismatch");
    const auto state=truth.at("states").doubles();const auto index=std::size_t(std::llround(p.time/f.sensorStep));
    if(index+1<f.basis.delays)throw std::runtime_error("insufficient prehistory");
    auto plant=core::createSystem(config);const auto sensor=plant->observableCatalog().at("far_field_acoustic_echo");
    const double initialEcho=sensor.evaluateReal(state.data()+2*index)[0];
    if(std::abs(initialEcho-rec.echo[index])>1e-12*std::max(1.,std::abs(initialEcho)))
        throw std::runtime_error("source sensor/configuration mismatch");
    std::vector<double> amplitudes{0,p.delta,-p.delta,p.delta/2,-p.delta/2,2*p.delta,-2*p.delta};
    std::vector<double> widths(amplitudes.size(),.01),steps(amplitudes.size(),f.sensorStep);
    if(p.extended)for(double width:{.005,.05})for(double sign:{1.,-1.}) {
        amplitudes.push_back(sign*p.delta);widths.push_back(width);steps.push_back(f.sensorStep);
    }
    if(p.refine)for(double amplitude:{0.,p.delta,-p.delta}){
        amplitudes.push_back(amplitude);widths.push_back(.01);steps.push_back(f.sensorStep/2);
    }
    std::vector<double> echoes,time;double replayError=0,replayEnergy=0;
    for(std::size_t branch=0;branch<amplitudes.size();++branch) {
        core::IntegrationOptions opts;opts.finalTime=14.;opts.outputStep=f.sensorStep;opts.maximumStep=steps[branch];
        opts.relativeTolerance=1e-9;opts.absoluteTolerance=1e-11;
        core::ScalarControlPolicy input;input.evaluation=core::ControlEvaluationMode::sampled;input.sampleTime=f.sensorStep;
        input.value=[&](std::size_t,double t,const double*){return t<widths[branch]-1e-10?amplitudes[branch]:0.;};
        const auto sol=core::integrateControlled(*plant,{state[2*index],state[2*index+1]},1,opts,input);
        if(branch==0)time=sol.time;
        if(sol.time!=time)throw std::runtime_error("branch timing mismatch");
        for(std::size_t k=0;k<sol.sampleCount;++k)echoes.push_back(sensor.evaluateReal(sol.state(0,k))[0]);
        for(std::size_t k=0;k<sol.sampleCount;++k) {
            const double expected=sol.time[k]<widths[branch]-1e-10?amplitudes[branch]:0.;
            if(std::abs(sol.controlValues[k]-expected)>1e-12)throw std::runtime_error("pulse timing error");
        }
        if(branch==0&&p.regime=="free")for(std::size_t k=0;k<sol.sampleCount&&index+k<rec.echo.size();++k) {
            replayError+=std::pow(echoes[k]-rec.echo[index+k],2);replayEnergy+=std::pow(rec.echo[index+k],2);
        }
    }
    Json meta=Json::Object{{"schema","dynsyslab.echo-pulse-acquisition/v1"},{"id",p.id},
        {"source",(folder/(p.source+".measurements.npz")).generic_string()},{"source_time",p.time},
        {"regime",p.regime},{"split",p.split},{"delta",p.delta},{"duration",14.},
        {"history_samples",double(f.basis.delays)},{"sensor_step",f.sensorStep},
        {"free_continuation_relative_error",relative(replayError,replayEnergy)},
        {"boundary","states used only to restart acquisition; only echo/input cross into analysis; intrinsic echo convention"}};
    core::NpzRecord a;
    a.arrays["time_since_pulse"]=core::NpyArray::fromDoubles({time.size()},time);
    a.arrays["echo_history"]=core::NpyArray::fromDoubles({f.basis.delays},std::vector<double>(rec.echo.begin()+index+1-f.basis.delays,rec.echo.begin()+index+1));
    a.arrays["echo_continuations"]=core::NpyArray::fromDoubles({amplitudes.size(),time.size()},echoes);
    a.arrays["input_amplitudes"]=core::NpyArray::fromDoubles({amplitudes.size()},amplitudes);
    a.arrays["pulse_widths"]=core::NpyArray::fromDoubles({widths.size()},widths);
    a.arrays["integration_max_steps"]=core::NpyArray::fromDoubles({steps.size()},steps);
    a.arrays["metadata_json"]=core::NpyArray::fromString(core::json::dump(meta));
    core::writeNpzAtomicNew(dest,a);return meta;
}
void acquireDense(const fs::path& root,const fs::path& out,const Frozen& f,const core::SystemConfiguration& config) {
    const auto folder=root/"archives/experiments/echo-input-regression-20261006-v1";
    auto plant=core::createSystem(config);const auto sensor=plant->observableCatalog().at("far_field_acoustic_echo");
    Json::Array rows;
    for(const std::string regime:{"mild","nonlinear"})for(int record=1;record<=8;++record) {
        const auto id="local-"+regime+"-"+std::to_string(record);
        const auto dest=out/("dense-"+id+".pulse.npz");
        if(fs::exists(dest)){rows.push_back(core::json::parse(core::readNpz(dest).metadataJson()));continue;}
        const auto measurements=readMeasurements(folder/(id+".measurements.npz"));
        const auto truth=core::readNpz(folder/(id+".truth.npz"));const auto state=truth.at("states").doubles();
        if(measurements.time!=truth.at("time").doubles())throw std::runtime_error("dense acquisition timing mismatch");
        const double delta=regime=="mild"?.01:.1;std::vector<double> indices,tails;
        const auto begin=std::size_t(std::llround(12/f.sensorStep)),end=std::size_t(std::llround((measurements.time.back()-6.02)/f.sensorStep));
        for(std::size_t k=begin;k+f.stride<=end;k+=4*f.stride) {
            indices.push_back(double(k));
            if(std::abs(sensor.evaluateReal(state.data()+2*k)[0]-measurements.echo[k])>1e-10)
                throw std::runtime_error("dense acquisition sensor mismatch");
            for(double sign:{1.,-1.}) {
                core::IntegrationOptions opts;opts.finalTime=.01;opts.outputStep=f.sensorStep;opts.maximumStep=f.sensorStep;
                opts.relativeTolerance=1e-9;opts.absoluteTolerance=1e-11;
                core::ScalarControlPolicy input;input.evaluation=core::ControlEvaluationMode::sampled;input.sampleTime=f.sensorStep;
                input.value=[&](std::size_t,double,const double*){return sign*delta;};
                const auto sol=core::integrateControlled(*plant,{state[2*k],state[2*k+1]},1,opts,input);
                if(sol.sampleCount!=f.stride+1)throw std::runtime_error("dense pulse sample count mismatch");
                for(std::size_t j=0;j<sol.sampleCount;++j)tails.push_back(sensor.evaluateReal(sol.state(0,j))[0]);
            }
        }
        Json meta=Json::Object{{"schema","dynsyslab.echo-dense-pulse/v1"},{"source",id},{"regime",regime},
            {"source_reference",(folder/(id+".measurements.npz")).generic_string()},
            {"split",record<=6?"train":"heldout"},{"delta",delta},{"points",double(indices.size())},
            {"stride_control_holds",4.},{"boundary","source echo histories referenced, no state crosses acquisition boundary"}};
        core::NpzRecord a;a.arrays["source_indices"]=core::NpyArray::fromDoubles({indices.size()},indices);
        a.arrays["echo_tails"]=core::NpyArray::fromDoubles({indices.size(),2,f.stride+1},tails);
        a.arrays["metadata_json"]=core::NpyArray::fromString(core::json::dump(meta));core::writeNpzAtomicNew(dest,a);rows.push_back(meta);
    }
    saveJson(out/"dense-acquisition.json",Json::Object{{"records",rows}});
}
}
