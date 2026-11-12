// Case-only client of the installed core. No fitting or controller implementation.
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/KoopmanValidation.hpp>
#include <dynsyslab/core/Json.hpp>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
namespace core=dynsyslab::core;
using J=core::json::Value;
using Path=std::filesystem::path;

void write(const Path& path,const J& value) {
    if(std::filesystem::exists(path)) throw std::runtime_error("Refusing to overwrite "+path.string());
    std::filesystem::create_directories(path.parent_path());
    std::ofstream out(path);out<<core::json::dump(value)<<'\n';
    if(!out) throw std::runtime_error("Unable to write "+path.string());
}
J complexValue(std::complex<double> z) {return J::Array{z.real(),z.imag()};}
J statistics(const core::KoopmanErrorStatistics& s) {
    return J::Object{{"pairs",double(s.samplePairs)},{"absolute_rms",s.absoluteRms},
        {"maximum_absolute",s.maximumAbsolute},{"relative_l2",s.relativeL2},
        {"median_pointwise_relative",s.medianPointwiseRelative},{"p95_pointwise_relative",s.p95PointwiseRelative}};
}
void prepare(const Path& specPath,const Path& out) {
    const auto spec=core::json::load(specPath.wstring());auto job=spec.at("solve");
    core::SystemConfiguration cfg;cfg.identifier=job.at("system").at("name").string();
    cfg.preset=job.at("system").at("preset").string();
    for(const auto& [key,value]:job.at("system").at("parameters").object())cfg.parameters[key]=value.number();
    const auto system=core::createSystem(cfg);
    const auto& input=spec.at("eigenfunction");
    std::vector<std::size_t> powers;std::vector<std::complex<double>> coefficients;
    for(const auto& row:input.at("exponents").array())for(const auto& value:row.array())powers.push_back(std::size_t(value.number()));
    for(const auto& row:input.at("coefficients_real_imag").array())coefficients.emplace_back(row.array().at(0).number(),row.array().at(1).number());
    const auto& rate=input.at("continuous_eigenvalue").array();
    auto polynomial=core::KoopmanSINDyResult::fromCoefficients(2,powers,coefficients,{{rate.at(0).number(),rate.at(1).number()}});
    polynomial.sampleTime=job.at("time").at("output_step").number();
    core::KoopmanModel model(polynomial,system->configuration());
    const auto modelPath=std::filesystem::absolute(out/"manual-eigenfunction.npz");
    if(std::filesystem::exists(modelPath))throw std::runtime_error("KROM already exists");
    std::filesystem::create_directories(out);
    core::savePolynomialKoopmanModel(modelPath,model,core::json::dump(J::Object{
        {"schema","dynsyslab.manual-polynomial-case/v1"},{"method","manual-polynomial"},
        {"label","RPE PDF manual eigenfunction"},{"data_role","krom"},{"source",spec.at("source")},
        {"integration_form","dimensionless"},{"delta_t",polynomial.sampleTime},
        {"training_data_policy","No refit; supplied coefficients kept unchanged"}}));
    const auto loaded=core::loadKoopmanModel(modelPath);
    double maximumRoundtrip=0;
    for(const std::vector<double> x: {std::vector<double>{1,0},{1.2,0},{1.3,0},{.9,-1},{1.1,1}})
        maximumRoundtrip=std::max(maximumRoundtrip,std::abs(model.evaluate(x.data())[0]-loaded.evaluate(x.data())[0]));
    if(maximumRoundtrip>1e-14)throw std::runtime_error("Archive changed supplied polynomial");
    job.at("controller").at("settings").at("krom")=modelPath.generic_string();
    write(out/"control-page2.json",job);
    auto refined=job;refined.at("label")="rpe-manual-pdf-page2-refined-20261022";
    refined.at("solver").at("max_step")=.0005;refined.at("solver").at("rtol")=1e-10;refined.at("solver").at("atol")=1e-12;
    write(out/"control-page2-refined.json",refined);
    auto free=job;free.at("label")="rpe-manual-pdf-free-validation-20261022";
    free.at("controller")=J::Object{{"kind","none"}};
    free.at("initial_conditions").at("values")=J::Array{J::Array{1.2,0.0}};
    free.at("time").at("stop")=5.0;free.at("archive").at("data_role")="training_data";
    write(out/"free-validation.json",free);
    J::Object groups;for(const auto& [key,value]:system->configuration().dimensionlessGroups)groups[key]=value;
    const double equilibrium[]{1,0},target[]{1.3,0};double drift[2];system->rhs(0,target,drift);
    const auto summary=J::Object{{"schema","dynsyslab.manual-rpe-preparation/v1"},
        {"krom",modelPath.generic_string()},{"roundtrip_max_error",maximumRoundtrip},
        {"continuous_eigenvalue",complexValue(model.continuousEigenvalues()[0])},
        {"discrete_eigenvalue",complexValue(model.discreteEigenvalues()[0])},
        {"psi_equilibrium",complexValue(model.evaluate(equilibrium)[0])},
        {"psi_target",complexValue(model.evaluate(target)[0])},
        {"target_static_input",-drift[1]/system->inputMatrix(target)[1]},
        {"dimensionless_groups",std::move(groups)}};
    write(out/"preparation.json",summary);std::cout<<core::json::dump(summary)<<'\n';
}
void analyze(const Path& modelPath,const Path& source,const Path& out,bool free) {
    const auto model=core::loadKoopmanModel(modelPath);const auto record=core::readNpz(source);
    const auto time=record.at("time").doubles(),values=record.at("values").doubles();
    const auto n=time.size();if(n<3||values.size()!=2*n)throw std::runtime_error("Expected one two-state trajectory");
    J::Object result{{"schema","dynsyslab.manual-rpe-analysis/v1"},{"source",source.generic_string()},{"samples",double(n)}};
    if(free) {
        const auto validation=core::validateKoopmanSemigroup(model,time,values,n,2,{.01,.1,1,5});
        const auto& v=validation.modes[0];result["one_step"]=statistics(v.oneStep);
        result["anchored_5_cycles"]=statistics(v.anchoredLongTerm);result["generator_finite_difference"]=statistics(v.generator);
        J::Array horizons;for(const auto& h:v.horizons)horizons.push_back(J::Object{{"horizon",h.horizon},{"statistics",statistics(h.statistics)}});
        result["horizons"]=horizons;
    } else {
        const auto control=record.at("control_values").doubles();
        if(control.size()!=n)throw std::runtime_error("Control sample mismatch");
        const double target[]{1.3,0};const auto desired=model.evaluate(target)[0];
        double radiusMin=values[0],radiusMax=values[0],peakVelocity=0,peakControl=0;
        double tailR=0,tailV=0,tailPsi=0,tailControl=0,tailCount=0;
        double entering1Percent=0,entering2Percent=0;
        for(std::size_t k=0;k<n;++k) {
            radiusMin=std::min(radiusMin,values[2*k]);radiusMax=std::max(radiusMax,values[2*k]);
            peakVelocity=std::max(peakVelocity,std::abs(values[2*k+1]));peakControl=std::max(peakControl,std::abs(control[k]));
            const double error=values[2*k]-target[0];
            if(std::abs(error)>.01*target[0])entering1Percent=k+1<n?time[k+1]:-1;
            if(std::abs(error)>.02*target[0])entering2Percent=k+1<n?time[k+1]:-1;
            if(time[k]>=30){tailR+=error*error;tailV+=values[2*k+1]*values[2*k+1];tailControl+=control[k];
                tailPsi+=std::norm(model.evaluate(&values[2*k])[0]-desired);++tailCount;}
        }
        result["radius_min"]=radiusMin;result["radius_max"]=radiusMax;result["peak_absolute_velocity"]=peakVelocity;
        result["final_state"]=J::Array{values[2*n-2],values[2*n-1]};result["peak_absolute_input"]=peakControl;
        result["initial_input"]=control.front();result["final_input"]=control.back();
        result["settling_time_1percent_of_target"]=entering1Percent;result["settling_time_2percent_of_target"]=entering2Percent;
        result["psi_target"]=complexValue(desired);result["psi_final"]=complexValue(model.evaluate(&values[2*n-2])[0]);
        if(tailCount){result["tail_30_40_radius_rmse"]=std::sqrt(tailR/tailCount);result["tail_30_40_velocity_rms"]=std::sqrt(tailV/tailCount);
            result["tail_30_40_psi_rmse"]=std::sqrt(tailPsi/tailCount);result["tail_30_40_mean_input"]=tailControl/tailCount;}
    }
    write(out,result);std::cout<<core::json::dump(result)<<'\n';
}
void compare(const Path& first,const Path& second,const Path& output) {
    const auto a=core::readNpz(first),b=core::readNpz(second);
    if(a.at("time").doubles()!=b.at("time").doubles())throw std::runtime_error("Comparison grids differ");
    const auto x=a.at("values").doubles(),y=b.at("values").doubles();
    const auto u=a.at("control_values").doubles(),v=b.at("control_values").doubles();
    if(x.size()!=y.size()||u.size()!=v.size())throw std::runtime_error("Comparison sizes differ");
    double radius=0,velocity=0,input=0;
    for(std::size_t i=0;i<x.size();i+=2){radius=std::max(radius,std::abs(x[i]-y[i]));velocity=std::max(velocity,std::abs(x[i+1]-y[i+1]));}
    for(std::size_t i=0;i<u.size();++i)input=std::max(input,std::abs(u[i]-v[i]));
    const J result=J::Object{{"schema","dynsyslab.manual-rpe-convergence/v1"},
        {"baseline",first.generic_string()},{"refined",second.generic_string()},
        {"maximum_absolute_radius_difference",radius},{"maximum_absolute_velocity_difference",velocity},
        {"maximum_absolute_input_difference",input}};
    write(output,result);std::cout<<core::json::dump(result)<<'\n';
}
int main(int argc,char** argv) {
    try {
        if(argc==4&&std::string(argv[1])=="prepare")prepare(argv[2],argv[3]);
        else if(argc==6&&std::string(argv[1])=="analyze")analyze(argv[2],argv[3],argv[4],std::string(argv[5])=="free");
        else if(argc==5&&std::string(argv[1])=="compare")compare(argv[2],argv[3],argv[4]);
        else throw std::runtime_error("prepare SPEC OUTPUT_DIR | analyze KROM SOLUTION OUTPUT_JSON free|controlled | compare BASELINE REFINED OUTPUT_JSON");
        return 0;
    } catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
