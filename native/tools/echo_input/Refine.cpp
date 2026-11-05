#include "Experiment.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

namespace echo_experiment {
namespace {
void saveLift(const std::filesystem::path& path,const Lifted& x) {
    if(std::filesystem::exists(path)) return;
    core::NpzRecord out;
    std::vector<double> xi;
    for(const auto& p:x.xi) xi.insert(xi.end(),p.begin(),p.end());
    out.arrays["time"]=core::NpyArray::fromDoubles({x.time.size()},x.time);
    out.arrays["xi"]=core::NpyArray::fromDoubles({x.time.size(),2},xi);
    out.arrays["psi"]=core::NpyArray::fromComplex({x.time.size()},x.psi);
    out.arrays["metadata_json"]=core::NpyArray::fromString(
        "{\"schema\":\"dynsyslab.echo-input-projection/v1\",\"basis\":\"frozen prior basis; no refit\"}");
    core::writeNpzAtomicNew(path,out);
}
Json fitScore(const core::HistoryInputMap& map,const std::vector<core::HistoryInputSample>& rows) {
    const auto s=core::scoreHistoryInputMap(map,rows);
    return Json::Object{{"relative_residual",s.relativeResidual},{"rmse",s.rmse},
        {"zero_input_baseline_rmse",s.unforcedBaselineRmse},{"pairs",double(rows.size())}};
}
Json domain(const Lifted& x,const std::vector<double>& lower,const std::vector<double>& upper) {
    double min0=1e300,max0=-1e300,min1=1e300,max1=-1e300,maxZ=0,q2=0,z2=0;
    std::size_t outside=0;
    for(std::size_t k=0;k<x.psi.size();++k) {
        const auto p=x.xi[k];
        min0=std::min(min0,p[0]);max0=std::max(max0,p[0]);
        min1=std::min(min1,p[1]);max1=std::max(max1,p[1]);
        outside+=(p[0]<lower[0] || p[0]>upper[0] || p[1]<lower[1] || p[1]>upper[1]);
        maxZ=std::max(maxZ,std::abs(x.psi[k])); z2+=std::norm(x.psi[k]);
    }
    return Json::Object{{"xi1_min",min0},{"xi1_max",max0},{"xi2_min",min1},{"xi2_max",max1},
        {"outside_training_box_fraction",double(outside)/x.psi.size()},
        {"psi_max",maxZ},{"psi_rms",std::sqrt(z2/x.psi.size())}};
}
}
void refineExperiment(const core::SystemConfiguration& config,const std::filesystem::path& folder,
    const Frozen& f,const core::NpzRecord& archive) {
    if(std::filesystem::exists(folder/"refinement.json"))
        throw std::invalid_argument("completed refinement is immutable");
    const auto train=archive.at("training_states").doubles();
    std::vector<double> lower(2,1e300),upper(2,-1e300);
    for(std::size_t k=0;k<train.size();++k) {
        lower[k%2]=std::min(lower[k%2],train[k]);upper[k%2]=std::max(upper[k%2],train[k]);
    }
    Json::Array datasets,fits,controls;
    for(const auto& name:std::vector<std::string>{"mild","nonlinear"}) {
        const double freq=name=="mild"?0.93:0.83;
        const auto target=readMeasurements(folder/("target-"+name+".measurements.npz"));
        const auto targetLift=lift(target,f);
        saveLift(folder/("target-"+name+".projection.npz"),targetLift);
        datasets.emplace_back(Json::Object{{"id","target-"+name},{"domain",domain(targetLift,lower,upper)}});
        std::vector<core::HistoryInputSample> local;
        const std::vector<double> amplitudes=name=="mild"?std::vector<double>{0.004,0.008,0.016,0.025}
            :std::vector<double>{0.18,0.24,0.36,0.42};
        int serial=0;
        const double stop=6+32/freq;
        for(double amp:amplitudes) for(double phase:{0.4,1.1}) {
            const auto id="local-"+name+"-"+std::to_string(++serial);
            if(!std::filesystem::exists(folder/(id+".measurements.npz"))) {
                std::cout<<"long calibration "<<id<<std::endl;
                std::size_t tick=0; double held=0;
                simulate(config,folder,id,f,stop+6,1.0,[&](double t,double) {
                    if(tick++%f.stride==0) held=t>=6-1e-9 && t<stop-1e-9
                        ?amp*std::sin(6.283185307179586*freq*(t-6)+phase):0;
                    return held;
                },Json::Object{{"role","frequency-local forced B identification; not eigenfunction training"},
                    {"amplitude_model",amp},{"frequency_cycles_per_tau",freq},{"phase",phase},
                    {"drive_start",6.0},{"drive_stop",stop}});
            }
            const auto x=lift(readMeasurements(folder/(id+".measurements.npz")),f);
            for(std::size_t k=0;k<x.pairs.size();++k)
                if(x.time[k]>=12 && x.time[k]<stop) local.push_back(x.pairs[k]);
            datasets.emplace_back(Json::Object{{"id",id},{"domain",domain(x,lower,upper)}});
        }
        for(std::size_t degree=0;degree<=1;++degree) {
            const auto id="local-"+name+"-degree-"+std::to_string(degree);
            const auto map=core::fitHistoryInputMap(local,f.multiplier,f.controlStep,degree,{1,1});
            const Json report=Json::Object{{"id",id},{"calibration",fitScore(map,local)},
                {"target_heldout_amplitude_phase",fitScore(map,targetLift.pairs)},
                {"design_condition",map.conditionNumber},{"rank",double(map.rank)}};
            fits.push_back(report); std::cout<<core::json::dump(report)<<std::endl;
            core::NpzRecord saved;
            saved.arrays["coefficients"]=core::NpyArray::fromComplex({map.coefficients.size()},map.coefficients);
            const auto canonical=map.canonicalRepresentation();
            if(const auto* c=std::get_if<core::ConstantControlMap>(&canonical)) {
                saved.arrays["matrix"]=core::NpyArray::fromComplex({1,1},c->matrix);
                saved.arrays["discrete_matrix"]=core::NpyArray::fromComplex({1,1},c->discreteMatrix);
                saved.arrays["sample_time"]=core::NpyArray::fromDoubles({1},{f.controlStep});
            } else {
                const auto& b=std::get<core::BilinearControlMap>(canonical);
                saved.arrays["operators"]=core::NpyArray::fromComplex({1,1,1},b.operators);
                saved.arrays["conjugate_operators"]=core::NpyArray::fromComplex({1,1,1},b.conjugateOperators);
                saved.arrays["offsets"]=core::NpyArray::fromComplex({1,1},b.offsets);
            }
            saved.arrays["metadata_json"]=core::NpyArray::fromString(core::json::dump(report));
            if(!std::filesystem::exists(folder/(id+".map.npz"))) core::writeNpzAtomicNew(folder/(id+".map.npz"),saved);
            for(double penalty:{100.0,1.0,0.01,0.0001}) {
                const auto attempt=id+"-klqr-r"+std::to_string(penalty);
                std::cout<<"control "<<attempt<<std::endl;
                try {
                    if(std::filesystem::exists(folder/(attempt+".measurements.npz"))) {
                        const auto actual=readMeasurements(folder/(attempt+".measurements.npz"));
                        controls.emplace_back(Json::Object{{"id",attempt},{"late",metrics(actual,target,30)}});
                    } else controls.push_back(controlTrial(config,folder,attempt,f,map,target,targetLift,penalty));
                    std::cout<<core::json::dump(controls.back())<<std::endl;
                } catch(const std::exception& e) { controls.emplace_back(Json::Object{{"id",attempt},{"failure",e.what()}}); }
            }
        }
    }
    // Per-waveform validation separates mild errors from large-polynomial outliers.
    Json::Array checks;
    for(int index=1;index<=9;++index) {
        const auto id="heldout-"+std::to_string(index);
        const auto x=lift(readMeasurements(folder/(id+".measurements.npz")),f);
        checks.emplace_back(Json::Object{{"id",id},{"domain",domain(x,lower,upper)}});
    }
    saveJson(folder/"refinement.json",Json::Object{{"calibration_datasets",datasets},
        {"maps",fits},{"controls",controls},{"broad_heldout_domains",checks},
        {"training_lower",Json::Array{lower[0],lower[1]}},{"training_upper",Json::Array{upper[0],upper[1]}}});
}
} // namespace echo_experiment
