// Run once before refactoring to freeze evidence, then compare after rebuilding.
// Uses public APIs only; it neither refits models nor changes scientific archives.
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/KoopmanValidation.hpp>
#include <dynsyslab/core/Json.hpp>
#include <fstream>
#include <iostream>
#include <cmath>
namespace c=dynsyslab::core;using J=c::json::Value;using P=std::filesystem::path;
int main(int argc,char** argv){try{
    if(argc!=4)throw std::runtime_error("ROOT BASELINE_JSON record|compare");
    const P root=argv[1],baseline=argv[2];std::vector<double> numbers;
    const auto model=c::loadKoopmanModel(root/"archives/experiments/rpe-forced-edmdc-validation-20261028-v1/all-modes.npz");
    auto append=[&](const auto& values){for(const auto z:values){numbers.push_back(std::real(z));numbers.push_back(std::imag(z));}};
    append(model.eigenfunctions().coefficients);append(model.continuousEigenvalues());append(model.discreteEigenvalues());
    for(int i=0;i<11;++i)for(int j=0;j<11;++j){const double x[]{.5+.2*i,-3.+.6*j};
        append(model.evaluate(x));append(model.exactControlMap(x));
        for(std::size_t m=0;m<model.liftedDimension();++m)append(c::gradientPolynomialEigenfunction(model.eigenfunctions(),m,x));}
    for(int i=0;i<5;++i){
        const auto archive=c::readNpz(root/("archives/solutions/rpe/rpe-viscous-echo-long-validation-five-ics-trajectory-"
            +std::to_string(i)+"__20261002T001403Z.npz"));
        const auto t=archive.at("time").doubles(),x=archive.at("values").doubles();std::vector<double> time,states;
        for(std::size_t k=0;k<t.size();k+=100){time.push_back(t[k]);states.insert(states.end(),{x[2*k],x[2*k+1]});}
        const auto check=c::validateKoopmanSemigroup(model,time,states,time.size(),2,{.1,1.,5.,10.});
        for(const auto& m:check.modes){numbers.insert(numbers.end(),{m.oneStep.relativeL2,m.anchoredLongTerm.relativeL2,
            m.generator.relativeL2,m.oneStepChangeRelativeL2});for(const auto& h:m.horizons)numbers.push_back(h.statistics.relativeL2);}
    }
    if(std::string(argv[3])=="record"){
        if(std::filesystem::exists(baseline))throw std::runtime_error("Baseline already exists");
        J::Array values;for(double n:numbers)values.emplace_back(n);std::ofstream out(baseline);
        out<<c::json::dump(J::Object{{"schema","dynsyslab.eigenfunction-api-parity/v1"},{"values",values}});
        if(!out)throw std::runtime_error("Baseline write failed");
    }else if(std::string(argv[3])=="compare"){
        const auto old=c::json::load(baseline.wstring()).at("values").array();
        if(old.size()!=numbers.size())throw std::runtime_error("Parity vector dimensions differ");
        for(std::size_t k=0;k<numbers.size();++k)if(!std::isfinite(numbers[k])||!std::isfinite(old[k].number()))
            throw std::runtime_error("Nonfinite parity quantity");
        double maximum=0;for(std::size_t k=0;k<numbers.size();++k)
            maximum=std::max(maximum,std::abs(numbers[k]-old[k].number())/std::max(1.,std::abs(old[k].number())));
        std::cout<<"Maximum scaled difference: "<<maximum<<'\n';
        if(maximum>1e-11)throw std::runtime_error("Numerical parity failed");
    }else throw std::runtime_error("Unknown mode");
    std::cout<<"Checked "<<numbers.size()<<" quantities across 15 modes and five trajectories\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
