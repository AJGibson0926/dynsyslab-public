// Same native joint B fitter, trained on paired causal sensitivities rather than
// contaminated total residuals. No state or equations enter this file.
#include "echo_pulse.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
namespace echo_pulse {
void fitAndControl(const fs::path& root,const fs::path& out,const std::string& window,
                   const Frozen& f,const core::KoopmanModel& bank,bool controls,bool dense) {
    std::map<std::string,ModeSamples> data;
    double cancellation=0;
    for(const auto& p:preparations()) {
        const auto a=core::readNpz(out/(window+"-"+p.id+".response.npz"));
        if(a.contains("states"))throw std::runtime_error("state-bearing regression input rejected");
        const auto z=a.at("psi").complexDoubles(),g=a.at("causal_response").complexDoubles();
        const auto xi=a.at("xi").doubles();
        const auto key=p.split+"-"+p.regime;if(data[key].empty())data[key].resize(3);
        for(std::size_t j=0;j<3;++j) {
            const auto mu=std::exp(bank.continuousEigenvalues()[j]*.01),next=mu*z[j]+g[3+j];
            cancellation=std::max(cancellation,std::abs((next-mu*z[j])-g[3+j])/std::max(1e-300,std::abs(g[3+j])));
            data[key][j].push_back({{xi[0],xi[1]},z[j],next,1.});
        }
    }
    if(cancellation>1e-5)throw std::runtime_error("paired sensitivity surrogate loses precision");
    if(dense)for(auto& item:denseSamples(root,out,window,f,bank))data[item.first]=std::move(item.second);
    data["train-pooled"].resize(3);
    for(const std::string regime:{"free","mild","nonlinear"})for(std::size_t j=0;j<3;++j) {
        const auto& source=data.at("train-"+regime)[j];
        auto& dest=data["train-pooled"][j];dest.insert(dest.end(),source.begin(),source.end());
    }
    Json::Array fits,trials;
    for(const std::string training:{"free","mild","nonlinear","pooled"})for(int mask=1;mask<8;++mask) {
        std::vector<std::size_t> selected;std::string tag;
        for(std::size_t j=0;j<3;++j)if(mask&(1<<j)){selected.push_back(j);tag+=std::to_string(j+1);}
        const auto fit=fitJointInput(data.at("train-"+training),selected,bank.continuousEigenvalues());
        auto model=bank.extract(selected);model.setControlMap("fitted",fit.map);
        Json::Object scores;
        for(const auto& item:data)scores[item.first]=scoreJointInput(fit,item.second,selected,bank.continuousEigenvalues());
        const auto name=window+(dense?"-dense-pulse-":"-pulse-")+training+"-m"+tag;
        Json result=Json::Object{{"window",window},{"training",training},{"subset",tag},
            {"condition",fit.condition},{"rank",double(fit.rank)},{"columns",double(fit.columns)},
            {"scores",scores},{"surrogate_max_relative_cancellation",cancellation},
            {"policy","G1 from paired pulses, each set jointly refitted; target histories score-only; frozen U and psi"}};
        fits.push_back(result);
        if(!fs::exists(out/(name+".fit.json")))saveJson(out/(name+".fit.json"),result);
        if(!fs::exists(out/(name+".krom.npz")))core::saveKoopmanModel(out/(name+".krom.npz"),model,
            core::json::dump(result),{{"sampled_joint_coefficients",core::NpyArray::fromComplex({selected.size(),fit.columns},fit.coefficients)}});
        if(!controls || (tag!="1"&&tag!="12") || (dense && training=="free"))continue;
        for(const std::string targetName:{"mild","nonlinear"}) {
            if(training!="free"&&training!="pooled"&&training!=targetName)continue;
            const auto inputs=root/"archives/experiments/echo-input-regression-20261006-v1";
            const auto target=readMeasurements(inputs/("target-"+targetName+".measurements.npz"));
            const auto reference=project(target,f);std::vector<double> scales(selected.size());
            for(std::size_t k=0;k<reference.time.size();++k){const auto z=model.evaluateCoordinates(reference.x.data()+2*k);
                for(std::size_t j=0;j<z.size();++j)scales[j]+=std::norm(z[j]);}
            for(auto& scale:scales)scale=std::sqrt(scale/reference.time.size());
            const double penalty=window=="short"?1.:100.;
            const auto id=name+"-"+targetName+"-r"+std::to_string(penalty);const auto dest=out/(id+".json");
            if(fs::exists(dest)){trials.push_back(readJson(dest));continue;}
            std::cout<<"CONTROL "<<id<<std::endl;
            Json trial;
            try {trial=runModelTrial(out,id,f,model,bank.systemConfiguration(),target,reference,scales,penalty);}
            catch(const std::exception& e){trial=Json::Object{{"id",id},{"failure",e.what()}};}
            trial.object()["training"]=training;trial.object()["subset"]=tag;trial.object()["target"]=targetName;
            saveJson(dest,trial);trials.push_back(trial);std::cout<<core::json::dump(trial)<<std::endl;
        }
    }
    const auto path=out/(window+(dense?"-dense":"")+(controls?"-pulse-controls.json":"-pulse-fits.json"));
    if(!fs::exists(path))saveJson(path,Json::Object{{"fits",fits},{"controls",trials},{"joint_fitter_tests",jointInputTests()}});
}
}
