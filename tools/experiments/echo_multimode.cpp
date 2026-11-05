#include "echo_multimode.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
namespace echo_multi {
Coordinates project(const Measurements& m,const Frozen& f,double begin,double end) {
    for(std::size_t k=1;k<m.time.size();++k)
        if(std::abs(m.time[k]-m.time[k-1]-f.sensorStep)>1e-8)
            throw std::invalid_argument("sensor spacing mismatch");
    const auto p=core::projectHistory(m.echo,f.basis);
    Coordinates out;
    const auto first=((f.basis.delays-1+f.stride-1)/f.stride)*f.stride;
    for(std::size_t k=first;k<m.time.size();k+=f.stride) {
        if(m.time[k]<begin-1e-9 || m.time[k]>end+1e-9)continue;
        const auto row=k-(f.basis.delays-1);
        out.time.push_back(m.time[k]);out.input.push_back(m.input[k]);
        out.x.push_back(p.coordinates[row*f.basis.rank]);
        out.x.push_back(p.coordinates[row*f.basis.rank+1]);
    }
    return out;
}
Json freeErrors(const core::KoopmanModel& model,const std::vector<Coordinates>& runs,
                double begin,double end) {
    const auto modes=model.liftedDimension();
    std::vector<double> pair(modes),one(modes),four(modes);
    for(const auto& r:runs) {
        Coordinates sub;
        for(std::size_t k=0;k<r.time.size();++k)if(r.time[k]>=begin-1e-9&&r.time[k]<=end+1e-9) {
            sub.time.push_back(r.time[k]);sub.x.push_back(r.x[2*k]);sub.x.push_back(r.x[2*k+1]);
        }
        const auto v=core::validateKoopmanSemigroup(model,sub.time,sub.x,sub.time.size(),2,{1.,4.});
        for(std::size_t j=0;j<modes;++j) {
            pair[j]+=v.modes[j].oneStep.relativeL2;
            one[j]+=v.modes[j].horizons[0].statistics.relativeL2;
            four[j]+=v.modes[j].horizons[1].statistics.relativeL2;
        }
    }
    Json::Array rows;
    for(std::size_t j=0;j<modes;++j)rows.emplace_back(Json::Object{
        {"mode",double(j)},{"pair",pair[j]/runs.size()},
        {"one_tau",one[j]/runs.size()},{"four_tau",four[j]/runs.size()}});
    return rows;
}
}
#ifndef DSL_ECHO_NO_MAIN
using namespace echo_multi;
int main(int argc,char** argv) {
    try {
        if(argc<4)throw std::invalid_argument("usage: echo_multimode ROOT NEW_OUTPUT long|short [--prepare|--joint]");
        const fs::path root=argv[1],out=argv[2]; const std::string window=argv[3];
        if(window=="--audit"){auditResults(root,out);return 0;}
        if(window!="long"&&window!="short")throw std::invalid_argument("window must be long or short");
        fs::create_directories(out);
        const bool joint=argc>=5 && std::string(argv[4])=="--joint";
        const auto inputs=root/"archives/experiments/echo-input-regression-20261006-v1";
        const auto refit=inputs/(window=="long"?"refit-xi2.krom.npz":"short-history-0p05/refit-xi2.krom.npz");
        const auto old=core::loadKoopmanModel(refit);
        const auto f=frozen(core::readNpz(refit),old.eigenfunctions());
        const auto config=old.systemConfiguration();
        const auto bankPath=(joint?root/"archives/experiments/echo-multimode-20261009-v2":out)/(window+"-bank.krom.npz");
        if(joint && !fs::exists(bankPath))throw std::runtime_error("joint comparison requires the exact prior bank");
        if(joint && !fs::exists(out/(window+"-joint-contract.json")))saveJson(out/(window+"-joint-contract.json"),
            Json::Object{{"schema","dynsyslab.echo-multimode-joint/v1"},{"bank_reference",bankPath.generic_string()},
                {"B_policy","refit all rows using all selected complex modes and conjugates; no fixed inherited rows"},
                {"ridge",1e-8},{"checks",jointInputTests()}});
        if(!fs::exists(bankPath)) {
            std::vector<Coordinates> training,validation;
            for(int k=1;k<=8;++k)training.push_back(project(readMeasurements(
                inputs/("free-train-"+std::to_string(k)+".measurements.npz")),f,35,80));
            for(int k=1;k<=6;++k)validation.push_back(project(readMeasurements(
                inputs/("free-validation-"+std::to_string(k)+".measurements.npz")),f,5,80));
            std::vector<std::vector<double>> runs;for(const auto& r:training)runs.push_back(r.x);
            core::PooledPolynomialEDMDOptions opts;opts.degree=5;opts.sampleTime=.01;opts.relativeRidge=1e-10;
            auto fit=core::fitPooledPolynomialEDMD(runs,2,opts);
            fit.eigenfunctions.coordinateProvider="causal-history";
            auto all=core::KoopmanModel(fit.eigenfunctions,config);
            // Only free-motion validation selects modes. Never select using control outcomes.
            const auto late=freeErrors(all,validation,35,80),early=freeErrors(all,validation,5,15);
            Json::Array spectrum;std::vector<std::size_t> selected;
            for(std::size_t j=0;j<all.liftedDimension();++j) {
                const auto l=all.continuousEigenvalues()[j];
                spectrum.emplace_back(Json::Object{{"mode",double(j)},{"lambda_re",l.real()},
                    {"lambda_im",l.imag()},{"late",late.array()[j]},{"early",early.array()[j]}});
            }
            for(const auto band:std::vector<std::pair<double,double>>{{4.,8.},{8.,16.},{16.,24.}}) {
                double best=1e300;std::size_t chosen=all.liftedDimension();
                for(std::size_t j=0;j<all.liftedDimension();++j) {
                    const auto l=all.continuousEigenvalues()[j];
                    if(l.real()>=-.00001 || l.real()<=-3 || l.imag()<band.first || l.imag()>=band.second)continue;
                    const auto& row=late.array()[j];const double s=row.at("one_tau").number()+row.at("four_tau").number();
                    if(s<best){best=s;chosen=j;}
                }
                if(chosen==all.liftedDimension())throw std::runtime_error("no stable candidate in predeclared band");
                selected.push_back(chosen);
            }
            Json::Array chosen;for(auto j:selected)chosen.push_back(spectrum[j]);
            const auto metadata=Json::Object{{"schema","dynsyslab.echo-multimode-bank/v1"},{"window",window},
                {"latent_coordinates",2.},{"history_samples",double(f.basis.delays)},
                {"history_span_tau",(f.basis.delays-1)*f.sensorStep},{"fit_start",35.},{"fit_end",80.},
                {"polynomial_degree",5.},{"spectrum",spectrum},{"selected",chosen},
                {"selection","stable positive-frequency modes; one per band [4,8), [8,16), [16,24); minimize free-validation 1+4 tau error on 35-80"},
                {"training","8 existing unforced training echoes; frozen U; no SVD; 6 reused free selection/validation records"},
                {"source_basis",refit.generic_string()}};
            core::saveKoopmanModel(out/(window+"-full-spectrum.krom.npz"),all,core::json::dump(metadata));
            core::saveKoopmanModel(bankPath,all.extract(selected),core::json::dump(metadata));
            saveJson(out/(window+"-bank.json"),metadata);
            std::cout<<"SELECTED "<<window<<" "<<core::json::dump(chosen)<<std::endl;
        }
        if(argc<5 || std::string(argv[4])!="--prepare")
            runControls(root,out,window,f,core::loadKoopmanModel(bankPath),config,joint);
        return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}
}
#endif
