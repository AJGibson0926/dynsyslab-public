// Case-only orchestration, diagnostics and archive preparation using native APIs.
// No new eigenfunction fitter, simulator, control law or error metric.
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/KoopmanValidation.hpp>
#include <dynsyslab/core/LQR.hpp>
#include <dynsyslab/core/Json.hpp>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
namespace c=dynsyslab::core;
using J=c::json::Value;
using Path=std::filesystem::path;
J complex(std::complex<double> z){return J::Array{z.real(),z.imag()};}
void write(const Path& p,const J& j){
    if(std::filesystem::exists(p))throw std::runtime_error("Refusing overwrite: "+p.string());
    std::filesystem::create_directories(p.parent_path());
    std::ofstream f(p);f<<c::json::dump(j)<<'\n';if(!f)throw std::runtime_error("write failed");
}
J accuracy(const c::KoopmanModel& m,const std::vector<double>& t,const std::vector<double>& x){
    const auto v=c::validateKoopmanSemigroup(m,t,x,t.size(),2,{.01,.1,1,5}).modes.at(0);
    return J::Object{{"one_step_relative_l2",v.oneStep.relativeL2},
        {"anchored_relative_l2",v.anchoredLongTerm.relativeL2},
        {"generator_fd_relative_l2",v.generator.relativeL2}};
}
// Map monomials by exponents, never assume the screenshot and SINDy use one order.
std::complex<double> align(const c::KoopmanSINDyResult& fitted,const c::KoopmanSINDyResult& manual){
    std::complex<double> numerator{};double denominator=0;
    for(std::size_t i=0;i<fitted.featureCount;++i){
        denominator+=std::norm(fitted.coefficients[i]);
        for(std::size_t j=0;j<manual.featureCount;++j)
            if(fitted.exponents[2*i]==manual.exponents[2*j]&&fitted.exponents[2*i+1]==manual.exponents[2*j+1])
                numerator+=std::conj(fitted.coefficients[i])*manual.coefficients[j];
    }
    return numerator/denominator;
}
J targetDiagnostic(const c::KoopmanModel& m){
    auto system=c::createSystem(m.systemConfiguration());const double x[]{1.3,0};double f[2];system->rhs(0,x,f);
    auto shared=std::make_shared<c::KoopmanModel>(m);
    c::KoopmanLQRController controller(shared,{2243,0,0,9.5},{.0001},
        [](double){return std::vector<double>{1.3,0};},"exact");
    const auto b=m.exactControlMap(x)[0],lambda=m.continuousEigenvalues()[0],psi=m.evaluate(x)[0];
    const auto gradient=c::gradientPolynomialEigenfunction(m.eigenfunctions(),0,x);
    const auto derivative=gradient[0]*f[0]+gradient[1]*f[1];
    return J::Object{{"psi",complex(psi)},{"u_from_normal_KLQR_at_target",controller.value(0,x)[0]},
        {"plant_free_acceleration_at_target",f[1]},
        {"required_static_input",-f[1]/system->inputMatrix(x)[1]},
        {"exact_Bpsi",complex(b)},{"generator_residual_at_target",complex(derivative-lambda*psi)},
        {"generator_relative_error_at_target",std::abs(derivative-lambda*psi)/std::abs(derivative)}};
}
int main(int argc,char** argv){try{
    if(argc!=5&&argc!=6)throw std::runtime_error("FIT TRAIN MANUAL_CONTROL_JOB OUTPUT_DIRECTORY [MODE_INDEX]");
    const Path output=std::filesystem::absolute(argv[4]);std::filesystem::create_directories(output);
    const auto all=c::loadKoopmanEigenfunctionModel(argv[1]);
    const auto train=c::readNpz(argv[2]);const auto time=train.at("time").doubles(),states=train.at("values").doubles();
    const auto job=c::json::load(Path(argv[3]).wstring());
    const auto manual=c::loadKoopmanModel(job.at("controller").at("settings").at("krom").string());
    std::vector<double> shortTime,shortStates;
    for(std::size_t k=0;k<time.size()&&time[k]<=5;++k){shortTime.push_back(time[k]);shortStates.push_back(states[2*k]);shortStates.push_back(states[2*k+1]);}
    std::size_t best=all.liftedDimension();double error=std::numeric_limits<double>::infinity();J::Array candidates,spectrum;
    for(std::size_t k=0;k<all.liftedDimension();++k){
        const auto lambda=all.continuousEigenvalues()[k];spectrum.push_back(complex(lambda));
        if(lambda.real()>=0||lambda.imag()<5||lambda.imag()>7)continue;
        const auto single=all.extract({k});const auto score=accuracy(single,time,states);
        candidates.push_back(J::Object{{"zero_based_index",double(k)},{"lambda",complex(lambda)},
            {"training_0_25",score},{"training_0_5",accuracy(single,shortTime,shortStates)}});
        if(score.at("anchored_relative_l2").number()<error){error=score.at("anchored_relative_l2").number();best=k;}
    }
    if(argc==6)best=std::stoul(argv[5]);
    if(best>=all.liftedDimension())throw std::runtime_error("No valid selected mode");
    const auto raw=all.extract({best});auto alignedFit=raw.eigenfunctions();const auto factor=align(alignedFit,manual.eigenfunctions());
    for(auto& v:alignedFit.coefficients)v*=factor;
    alignedFit.trainingValues.clear();alignedFit.residualRms.clear();alignedFit.residualMaximum.clear();
    const c::KoopmanModel aligned(alignedFit,raw.systemConfiguration());
    J::Array runs;
    for(const auto& name:{std::string("raw"),std::string("gauge-aligned")}){
        const auto& model=name=="raw"?raw:aligned;const auto krom=output/(name+".npz");
        if(std::filesystem::exists(krom))throw std::runtime_error("KROM exists");
        c::savePolynomialKoopmanModel(krom,model,c::json::dump(J::Object{
            {"method","case-native-extracted-sindy"},{"source_krom",argv[1]},
            {"source_training",argv[2]},{"selected_mode",double(best)},{"label",name},
            {"gauge_multiplier",complex(name=="raw"?std::complex<double>{1,0}:factor)}}));
        auto controlled=job;controlled.at("label")="rpe-pdf-retrained-mode"+std::to_string(best)+"-"+name+"-20261022";
        controlled.at("controller").at("settings").at("krom")=krom.generic_string();
        write(output/("control-"+name+".json"),controlled);
        runs.push_back(J::Object{{"name",name},{"krom",krom.generic_string()},
            {"free_0_5",accuracy(model,shortTime,shortStates)},
            {"free_0_25",accuracy(model,time,states)},{"target_diagnostic",targetDiagnostic(model)}});
    }
    double minR=states[0],maxR=states[0];for(std::size_t k=0;k<time.size();++k){minR=std::min(minR,states[2*k]);maxR=std::max(maxR,states[2*k]);}
    const J result=J::Object{{"selection_policy",argc==6?"Explicit mode cross-check; no control outcome used":"minimum anchored 0..25 free-training error among stable modes with 5<Im(lambda)<7; no control outcome used"},
        {"training",argv[2]},{"fit",argv[1]},{"all_spectrum",spectrum},{"candidates",candidates},
        {"selected_zero_based_index",double(best)},{"selected_lambda",complex(raw.continuousEigenvalues()[0])},
        {"gauge_multiplier",complex(factor)},{"training_radius_range",J::Array{minR,maxR}},
        {"manual_free_0_5",accuracy(manual,shortTime,shortStates)},
        {"manual_free_0_25",accuracy(manual,time,states)},
        {"manual_target_diagnostic",targetDiagnostic(manual)},{"runs",runs}};
    write(output/"preparation.json",result);std::cout<<c::json::dump(result)<<'\n';return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
