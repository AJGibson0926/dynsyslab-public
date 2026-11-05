// Isolated audit executable. Links the existing native DLLs and Experiment.cpp.
// No model fitting or controller changes. Only acquireProbes may access a plant.
#include "Experiment.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>

using namespace echo_experiment;
namespace fs = std::filesystem;
namespace {
double ratio(double a,double b) { return std::sqrt(a/std::max(b,1e-300)); }
core::HistoryInputMap loadMap(const fs::path& path,const Frozen& f) {
    core::HistoryInputMap map;
    map.coefficients=core::readNpz(path).at("coefficients").complexDoubles();
    map.degree=map.coefficients.size()==1?0:1;
    map.sampleTime=f.controlStep; map.multiplier=f.multiplier;
    return map;
}
struct Case { std::string id; fs::path folder; Frozen f; };
struct Probe {
    double time{},delta{},knownInput{};
    std::vector<double> history,inputs;
    std::vector<std::vector<double>> tails;
};
// Acquisition-only simulation adapter. State never crosses this boundary.
// All probes share a complete measured prehistory and branch for just one hold.
std::vector<Probe> acquireProbes(const core::SystemConfiguration& config,
    const Measurements& recording,const Frozen& f,double frequency,double delta,
    const fs::path& path) {
    auto plant=core::createSystem(config);
    const auto sensor=plant->observableCatalog().at("far_field_acoustic_echo");
    core::IntegrationOptions opts;
    opts.finalTime=32; opts.outputStep=f.sensorStep; opts.maximumStep=f.sensorStep;
    opts.relativeTolerance=1e-9; opts.absoluteTolerance=1e-11;
    core::ScalarControlPolicy replay;
    replay.evaluation=core::ControlEvaluationMode::sampled; replay.sampleTime=f.controlStep;
    replay.value=[&](std::size_t,double t,const double*) {
        return recording.input.at(std::llround(t/f.sensorStep));
    };
    const auto prefix=core::integrateControlled(*plant,{1,0},1,opts,replay);
    std::vector<double> echo(prefix.sampleCount);
    double replayError=0,replayEnergy=0;
    for(std::size_t k=0;k<echo.size();++k) {
        echo[k]=sensor.evaluateReal(prefix.state(0,k))[0];
        replayError+=std::pow(echo[k]-recording.echo[k],2); replayEnergy+=std::pow(recording.echo[k],2);
    }
    std::vector<Probe> probes;
    std::vector<double> histories,tails,inputs,times;
    for(int j=0;j<12;++j) {
        const auto tick=static_cast<std::size_t>(std::llround((30+j/(12*frequency))/f.controlStep));
        const auto index=tick*f.stride;
        Probe p; p.time=prefix.time[index]; p.delta=delta; p.knownInput=recording.input[index];
        p.history.assign(echo.begin()+index+1-f.basis.delays,echo.begin()+index+1);
        p.inputs={0,delta,-delta,2*delta,-2*delta,delta/2,-delta/2,p.knownInput};
        const auto state=prefix.state(0,index);
        for(double input:p.inputs) {
            core::ScalarControlPolicy held;
            held.evaluation=core::ControlEvaluationMode::sampled; held.sampleTime=f.sensorStep;
            held.value=[input](std::size_t,double,const double*) { return input; };
            auto nextOpts=opts; nextOpts.initialTime=p.time; nextOpts.finalTime=p.time+f.controlStep;
            const auto next=core::integrateControlled(*plant,{state[0],state[1]},1,nextOpts,held);
            std::vector<double> tail;
            for(std::size_t k=1;k<=f.stride;++k) tail.push_back(sensor.evaluateReal(next.state(0,k))[0]);
            p.tails.push_back(tail); tails.insert(tails.end(),tail.begin(),tail.end());
        }
        histories.insert(histories.end(),p.history.begin(),p.history.end());
        inputs.insert(inputs.end(),p.inputs.begin(),p.inputs.end()); times.push_back(p.time);
        probes.push_back(std::move(p));
    }
    core::NpzRecord out;
    out.arrays["time"]=core::NpyArray::fromDoubles({12},times);
    out.arrays["echo_history"]=core::NpyArray::fromDoubles({12,f.basis.delays},histories);
    out.arrays["input"]=core::NpyArray::fromDoubles({12,8},inputs);
    out.arrays["echo_tail"]=core::NpyArray::fromDoubles({12,8,f.stride},tails);
    out.arrays["metadata_json"]=core::NpyArray::fromString(core::json::dump(Json::Object{
        {"role","matched-history causal diagnostic, not training or independent validation"},
        {"replay_echo_relative_l2",ratio(replayError,replayEnergy)},
        {"boundary","only acquisition knows the simulator; diagnostic inputs contain echo, time and input only"}}));
    core::writeNpzAtomicNew(path,out);
    std::cout<<"probe acquisition replay error "<<ratio(replayError,replayEnergy)<<std::endl;
    return probes;
}
Complex psiHistory(const std::vector<double>& history,const Frozen& f) {
    core::CausalHistoryProjector projector(f.basis);
    for(double e:history) projector.append(e);
    const auto xi=projector.projectLatest();
    return core::evaluatePolynomialEigenfunction(f.polynomial,0,xi->data());
}
Json matched(const std::vector<Probe>& probes,const Case& c,const core::HistoryInputMap& map,
    const fs::path& output) {
    double drift=0,actuation=0,total=0,bErr=0,bEnergy=0,halfErr=0,doubleErr=0,affineErr=0,uError=0,uEnergy=0;
    std::vector<double> times,known,inverted,drifts,acts,bMag,mapMag;
    std::vector<Complex> bs,qs;
    std::vector<core::HistoryInputSample> train,check;
    std::size_t phase=0;
    for(const auto& p:probes) {
        const auto z=psiHistory(p.history,c.f);
        std::vector<Complex> after;
        for(const auto& tail:p.tails) {
            auto h=p.history; h.insert(h.end(),tail.begin(),tail.end());
            after.push_back(psiHistory(h,c.f));
        }
        const auto q0=after[0]-c.f.multiplier*z;
        const auto b=(after[1]-after[2])/(2*p.delta);
        const auto half=(after[5]-after[6])/p.delta;
        const auto twice=(after[3]-after[4])/(4*p.delta);
        const auto bu=after[7]-after[0], q=after[7]-c.f.multiplier*z;
        drift+=std::norm(q0); actuation+=std::norm(bu); total+=std::norm(q);
        bErr+=std::norm(map.evaluate(z)-b); bEnergy+=std::norm(b);
        halfErr+=std::norm(half-b); doubleErr+=std::norm(twice-b);
        affineErr+=std::norm(bu-b*p.knownInput);
        times.push_back(p.time); known.push_back(p.knownInput);
        inverted.push_back((std::conj(b)*(q-q0)).real()/std::max(std::norm(b),1e-300));
        uError+=std::pow(inverted.back()-p.knownInput,2); uEnergy+=p.knownInput*p.knownInput;
        drifts.push_back(std::abs(q0)); acts.push_back(std::abs(bu));
        bMag.push_back(std::abs(b)); mapMag.push_back(std::abs(map.evaluate(z)));
        bs.push_back(b); qs.push_back(q0);
        // Same bilinear fitter, but paired differences remove autonomous drift.
        // Every third phase is held out. This is a tiny local capacity check.
        (phase++%3==0?check:train).push_back({{},z,c.f.multiplier*z+b,1.});
    }
    core::NpzRecord out;
    for(const auto& kv:std::vector<std::pair<std::string,std::vector<double>>>{
        {"time",times},{"known_input",known},{"paired_inverse",inverted},{"drift",drifts},
        {"actuation",acts},{"causal_B_magnitude",bMag},{"fitted_B_magnitude",mapMag}})
        out.arrays[kv.first]=core::NpyArray::fromDoubles({kv.second.size()},kv.second);
    out.arrays["causal_B"]=core::NpyArray::fromComplex({bs.size()},bs);
    out.arrays["zero_input_residual"]=core::NpyArray::fromComplex({qs.size()},qs);
    core::writeNpzAtomicNew(output,out);
    const auto pairedMap=core::fitHistoryInputMap(train,c.f.multiplier,c.f.controlStep,1,{1,1});
    const auto fit=core::scoreHistoryInputMap(pairedMap,train),test=core::scoreHistoryInputMap(pairedMap,check);
    return Json::Object{{"zero_input_drift_over_actual_input_effect",ratio(drift,actuation)},
        {"drift_over_total_residual",ratio(drift,total)},
        {"fitted_B_relative_error_to_causal_B",ratio(bErr,bEnergy)},
        {"causal_B_half_delta_change",ratio(halfErr,bEnergy)},
        {"causal_B_double_delta_change",ratio(doubleErr,bEnergy)},
        {"paired_inverse_driver_relative_error",ratio(uError,uEnergy)},
        {"paired_bilinear_training_relative_error",fit.relativeResidual},
        {"paired_bilinear_heldout_phase_relative_error",test.relativeResidual},
        {"paired_bilinear_rank",double(pairedMap.rank)},{"paired_bilinear_condition",pairedMap.conditionNumber},
        {"local_input_affine_residual",ratio(affineErr,actuation)}};
}
std::vector<Probe> loadProbes(const fs::path& path,double delta) {
    const auto a=core::readNpz(path);
    if(a.contains("states")) throw std::invalid_argument("diagnostic rejects state-bearing probes");
    const auto time=a.at("time").doubles(),history=a.at("echo_history").doubles();
    const auto input=a.at("input").doubles(),tails=a.at("echo_tail").doubles();
    const auto delays=history.size()/time.size(),stride=tails.size()/(time.size()*8);
    std::vector<Probe> probes;
    for(std::size_t k=0;k<time.size();++k) {
        Probe p; p.time=time[k]; p.delta=delta; p.knownInput=input[k*8+7];
        p.history.assign(history.begin()+k*delays,history.begin()+(k+1)*delays);
        p.inputs.assign(input.begin()+k*8,input.begin()+(k+1)*8);
        for(std::size_t j=0;j<8;++j) p.tails.emplace_back(tails.begin()+(k*8+j)*stride,tails.begin()+(k*8+j+1)*stride);
        probes.push_back(std::move(p));
    }
    return probes;
}
Json releaseCheck(const Measurements& m,const Case& c,const fs::path& output) {
    const auto x=lift(m,c.f); Json::Array metrics;
    std::vector<double> time,residual,magnitude;
    for(std::size_t k=0;k<x.pairs.size();++k) if(x.time[k]>=30) {
        time.push_back(x.time[k]-30); magnitude.push_back(std::abs(x.psi[k]));
        residual.push_back(std::abs(x.psi[k+1]-c.f.multiplier*x.psi[k]));
    }
    for(auto band:std::vector<std::pair<double,double>>{{0,.05},{.05,.5},{.5,1},{1,4.3215},{4.3215,8},{8,12}}) {
        double err=0,energy=0;
        for(std::size_t k=0;k<time.size();++k) if(time[k]>=band.first&&time[k]<band.second) {
            err+=residual[k]*residual[k]; energy+=magnitude[k]*magnitude[k];
        }
        metrics.emplace_back(Json::Object{{"time_since_release_from",band.first},
            {"time_since_release_to",band.second},{"pair_relative_error",ratio(err,energy)}});
    }
    core::NpzRecord out;
    out.arrays["time_since_release"]=core::NpyArray::fromDoubles({time.size()},time);
    out.arrays["residual_magnitude"]=core::NpyArray::fromDoubles({time.size()},residual);
    out.arrays["psi_magnitude"]=core::NpyArray::fromDoubles({time.size()},magnitude);
    core::writeNpzAtomicNew(output,out); return metrics;
}
Json inverseTarget(const Measurements& m,const Lifted& x,const Frozen& f,
    const core::HistoryInputMap& map,const fs::path& output) {
    double q2=0,perp2=0,pred2=0,u2=0,uerr2=0,urec2=0,dot=0,z2=0;
    double minB=1e300,maxB=0; std::size_t saturated=0,n=0;
    std::vector<double> time,known,inverse,perp,magq,magb;
    for(std::size_t k=0;k<x.pairs.size();++k) {
        const auto& s=x.pairs[k]; const auto q=s.nextPsi-f.multiplier*s.psi,b=map.evaluate(s.psi);
        const double u=(std::conj(b)*q).real()/std::max(std::norm(b),1e-300);
        if(x.time[k]>=30) {
            q2+=std::norm(q); perp2+=std::norm(q-b*u); pred2+=std::norm(q-b*s.input);
            u2+=s.input*s.input; uerr2+=std::pow(u-s.input,2); urec2+=u*u; dot+=u*s.input;
            z2+=std::norm(s.psi); minB=std::min(minB,std::abs(b)); maxB=std::max(maxB,std::abs(b));
            saturated+=std::abs(u)>.6; ++n;
        }
        time.push_back(x.time[k]); known.push_back(s.input); inverse.push_back(u);
        perp.push_back(std::abs(q-b*u)); magq.push_back(std::abs(q)); magb.push_back(std::abs(b));
    }
    core::NpzRecord out;
    for(const auto& kv:std::vector<std::pair<std::string,std::vector<double>>>{
        {"time",time},{"known_input",known},{"inverse_input",inverse},{"perpendicular_residual",perp},
        {"total_residual",magq},{"B_magnitude",magb}})
        out.arrays[kv.first]=core::NpyArray::fromDoubles({kv.second.size()},kv.second);
    core::writeNpzAtomicNew(output,out);
    return Json::Object{{"late_inverse_driver_relative_error",ratio(uerr2,u2)},
        {"late_inverse_driver_cosine",dot/std::sqrt(std::max(u2*urec2,1e-300))},
        {"late_perpendicular_residual_fraction",ratio(perp2,q2)},
        {"late_known_driver_model_residual",ratio(pred2,q2)},
        {"late_residual_over_psi",ratio(q2,z2)},{"inverse_outside_actuator_fraction",double(saturated)/n},
        {"minimum_B",minB},{"maximum_B",maxB}};
}
Json freeValidation(const Case& c,const fs::path& data) {
    Json::Array intervals;
    for(auto interval:std::vector<std::pair<double,double>>{{5,15},{15,35},{50,75}}) {
        double residual=0,energy=0,one=0,oneEnergy=0,four=0,fourEnergy=0;
        for(int j=1;j<=6;++j) {
            const auto x=lift(readMeasurements(data/("free-validation-"+std::to_string(j)+".measurements.npz")),c.f);
            for(std::size_t k=0;k+400<x.psi.size();++k) if(x.time[k]>=interval.first&&x.time[k]<interval.second) {
                residual+=std::norm(x.psi[k+1]-c.f.multiplier*x.psi[k]); energy+=std::norm(x.psi[k]);
                one+=std::norm(x.psi[k+100]-std::pow(c.f.multiplier,100)*x.psi[k]); oneEnergy+=std::norm(x.psi[k+100]);
                four+=std::norm(x.psi[k+400]-std::pow(c.f.multiplier,400)*x.psi[k]); fourEnergy+=std::norm(x.psi[k+400]);
            }
        }
        intervals.emplace_back(Json::Object{{"from",interval.first},{"to",interval.second},
            {"pair_relative_error",ratio(residual,energy)},{"one_tau_relative_error",ratio(one,oneEnergy)},
            {"four_tau_relative_error",ratio(four,fourEnergy)}});
    }
    return intervals;
}
Json parity(const Measurements& m,const Frozen& f) {
    const auto all=core::projectHistory(m.echo,f.basis);
    core::CausalHistoryProjector online(f.basis);
    double err=0,energy=0,maxerr=0; std::size_t checks=0;
    for(std::size_t k=0;k<m.echo.size();++k) {
        online.append(m.echo[k]);
        if(k+1<f.basis.delays||k%1000) continue;
        const auto x=online.projectLatest(); const auto row=k+1-f.basis.delays;
        for(std::size_t j=0;j<f.basis.rank;++j) {
            const auto a=all.coordinates[row*f.basis.rank+j],d=(*x)[j]-a;
            err+=d*d; energy+=a*a; maxerr=std::max(maxerr,std::abs(d));
        }
        ++checks;
    }
    return Json::Object{{"relative_l2",ratio(err,energy)},{"maximum_absolute",maxerr},{"windows_checked",double(checks)}};
}
} // namespace
int main(int argc,char**argv) {
    try {
        if(argc!=3&&argc!=4) throw std::invalid_argument("usage: audit ROOT OUTPUT [--skip-controls|--paired-only]");
        const bool skipControls=argc==4&&std::string(argv[3])=="--skip-controls";
        const bool pairedOnly=argc==4&&std::string(argv[3])=="--paired-only";
        if(argc==4&&!skipControls&&!pairedOnly) throw std::invalid_argument("unknown option");
        const fs::path root=argv[1],output=argv[2],data=root/"archives/experiments/echo-input-regression-20261006-v1";
        fs::create_directories(output);
        if(fs::exists(output/"diagnosis.json")) throw std::runtime_error("completed diagnostic is immutable");
        std::vector<Case> cases;
        for(const auto& window:std::vector<std::string>{"long","short"}) for(int d:{2,4}) {
            const auto folder=window=="long"?data:data/"short-history-0p05";
            const auto path=folder/("refit-xi"+std::to_string(d)+".krom.npz");
            const auto model=core::loadKoopmanModel(path);
            cases.push_back({window+"-xi"+std::to_string(d),folder,frozen(core::readNpz(path),model.eigenfunctions())});
        }
        const auto config=core::loadKoopmanModel(data/"refit-xi2.krom.npz").systemConfiguration();
        if(pairedOnly) {
            Json::Array rows;
            for(const auto& target:std::vector<std::string>{"mild","nonlinear"}) {
                const auto probes=loadProbes(root/"archives/experiments/echo-control-diagnosis-20261007-v2"/(target+"-probes.npz"),target=="mild"?.01:.1);
                for(const auto& c:cases) {
                    const auto name="refit-xi"+std::to_string(c.f.polynomial.stateDimension)+"-"+target+"-local";
                    const auto map=loadMap(c.folder/(name+".map.npz"),c.f);
                    rows.emplace_back(Json::Object{{"case",c.id},{"target",target},
                        {"metrics",matched(probes,c,map,output/(c.id+"-"+target+"-paired.npz"))}});
                }
            }
            saveJson(output/"diagnosis.json",Json::Object{{"paired_capacity",rows}}); return 0;
        }
        Json::Array results,free,parities,resumed,releases;
        for(const auto& c:cases) {
            std::cout<<"free validation "<<c.id<<std::endl;
            free.emplace_back(Json::Object{{"case",c.id},{"intervals",freeValidation(c,data)}});
        }
        for(const auto& target:std::vector<std::string>{"mild","nonlinear"}) {
            const auto m=readMeasurements(data/("target-"+target+".measurements.npz"));
            const auto releaseId=target+"-release";
            simulate(config,output,releaseId,cases[0].f,42,1,[&](double t,double) {
                return t<30-1e-9?m.input.at(std::llround(t/cases[0].f.sensorStep)):0.;
            },Json::Object{{"role","diagnostic: replay target then release to zero input at tau=30"}});
            const auto released=readMeasurements(output/(releaseId+".measurements.npz"));
            const auto probes=acquireProbes(config,m,cases[0].f,target=="mild"?.93:.83,
                target=="mild"?.01:.1,output/(target+"-probes.npz"));
            for(const auto& c:cases) {
                std::cout<<"target "<<c.id<<" "<<target<<std::endl;
                const auto x=lift(m,c.f);
                releases.emplace_back(Json::Object{{"case",c.id},{"target",target},
                    {"intervals",releaseCheck(released,c,output/(c.id+"-"+target+"-release.npz"))}});
                parities.emplace_back(Json::Object{{"case",c.id},{"target",target},{"check",parity(m,c.f)}});
                for(const auto& family:std::vector<std::string>{"local","broad"}) {
                    const auto name="refit-xi"+std::to_string(c.f.polynomial.stateDimension)+"-"+target+"-"+family;
                    const auto map=loadMap(c.folder/(name+".map.npz"),c.f);
                    const auto id=c.id+"-"+target+"-"+family;
                    const auto inv=inverseTarget(m,x,c.f,map,output/(id+"-inverse.npz"));
                    const auto matchedResult=matched(probes,c,map,output/(id+"-paired.npz"));
                    results.emplace_back(Json::Object{{"id",id},{"inverse",inv},{"matched_history",matchedResult}});
                    std::cout<<core::json::dump(results.back())<<std::endl;
                    if(c.id=="short-xi4"&&target=="nonlinear"&&family=="local"&&!skipControls) {
                        for(double penalty:{100.,10.,1.,.01}) {
                            const auto trial=name+"-klqr-r"+std::to_string(penalty);
                            std::cout<<"resume storage-interrupted "<<trial<<std::endl;
                            try { resumed.push_back(controlTrial(config,output,trial,c.f,map,m,x,penalty)); }
                            catch(const std::exception& e) { resumed.emplace_back(Json::Object{{"id",trial},{"failure",e.what()}}); }
                            std::cout<<core::json::dump(resumed.back())<<std::endl;
                        }
                    }
                }
            }
        }
        saveJson(output/"diagnosis.json",Json::Object{{"results",results},{"free_validation",free},
            {"projection_parity",parities},{"resumed_controls",resumed},{"release_checks",releases},
            {"contract","existing DLL computations, frozen U/psi/B/KLQR; no core changes or refitting"}});
        return 0;
    } catch(const std::exception& e) {std::cerr<<e.what()<<std::endl;return 1;}
}
