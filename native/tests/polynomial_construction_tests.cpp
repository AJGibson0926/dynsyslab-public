#include <dynsyslab/core/KoopmanArchive.hpp>

#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
namespace core=dynsyslab::core;
void require(bool ok,const char* message) {if(!ok) throw std::runtime_error(message);}
template<class F> void rejects(F f) {
    bool failed=false;try{f();}catch(const std::invalid_argument&){failed=true;}
    require(failed,"invalid manual polynomial accepted");
}
}
int main() {
    std::filesystem::path path;
    try {
        const std::vector<std::size_t> powers{0,0,1,0,0,1};
        const std::vector<std::complex<double>> coefficients{{2,0},{1,0},{0,1}};
        const std::complex<double> lambda(-.2,3);
        const auto rates=[](std::vector<std::complex<double>> values){
            return core::Spectrum::fromContinuousRates(std::move(values));};
        auto polynomial=core::polynomialEigenfunctions(2,powers,coefficients,rates({lambda}));
        const double q[]{3,4};
        require(std::abs(polynomial.functions.evaluate(0,q)-std::complex<double>(5,4))<1e-14,
            "constant/complex coefficients not preserved");
        require(polynomial.diagnostics.sampleCount==0&&polynomial.diagnostics.directFitCount==0
            &&!polynomial.spectrum.hasResiduals(),"invented fit evidence");
        const auto gradient=polynomial.functions.gradient(0,q);
        require(gradient==std::vector<std::complex<double>>({{1,0},{0,1}}),"manual analytic gradient failed");
        rejects([&]{(void)polynomial.spectrum.discreteMultipliers();});
        polynomial.spectrum=core::Spectrum::fromContinuousRates({lambda},.1);
        const auto mu=polynomial.spectrum.discreteMultipliers()[0];
        require(std::abs(mu-std::exp(lambda*.1))<1e-14,"continuous to discrete failed");
        require(std::abs(polynomial.spectrum.discreteMultipliers(.2)[0]-mu*mu)<1e-14,"new dt conversion failed");
        auto historyDomain=core::FunctionDomain::supplied(2);
        historyDomain.providerIdentity="causal-history";
        historyDomain.names={"xi_1","xi_2"};
        auto sampled=core::polynomialEigenfunctions(2,powers,coefficients,
            core::Spectrum::fromDiscreteMultipliers({mu},.1),{},{},historyDomain);
        require(std::abs(sampled.spectrum.continuousRates()[0]-lambda)<1e-14,"discrete to continuous failed");
        auto aliased=core::polynomialEigenfunctions(2,powers,coefficients,
            core::Spectrum::fromDiscreteMultipliers({mu},.1,{1}));
        require(std::abs(aliased.spectrum.discreteMultipliers()[0]-mu)<1e-13,"branch did not preserve multiplier");
        require(std::abs(aliased.spectrum.continuousRates()[0].imag()-3-20*std::acos(-1.0))<1e-12,"branch rate incorrect");
        rejects([&]{(void)core::Spectrum::fromDiscreteMultipliers({0},.1).continuousRates();});
        rejects([&]{(void)core::polynomialEigenfunctions(2,powers,coefficients,
            core::Spectrum::fromDiscreteMultipliers({mu},0));});
        rejects([&]{(void)core::polynomialEigenfunctions(2,powers,{1},rates({lambda}));});
        rejects([&]{(void)core::polynomialEigenfunctions(2,powers,{0,0,0},rates({lambda}));});
        rejects([&]{(void)core::polynomialEigenfunctions(2,{0,0,0,0},{1,2},rates({lambda}));});
        rejects([&]{(void)core::polynomialEigenfunctions(2,powers,coefficients,rates({lambda}),{0,0},{1,0});});
        rejects([&]{(void)core::polynomialEigenfunctions(2,powers,coefficients,
            rates({{std::numeric_limits<double>::quiet_NaN(),0}}));});
        core::SystemConfiguration system; system.identifier="duffing";
        const auto direct=core::KoopmanModel::fromPolynomialCoefficients(system,2,powers,coefficients,{lambda});
        require(std::abs(direct.evaluate(q)[0]-std::complex<double>(5,4))<1e-14,"manual KROM evaluation differs");
        core::KoopmanModel model(sampled,system);
        require(model.requiresExternalCoordinates(),"history coordinates became physical state");
        require(std::abs(model.evaluateCoordinates(q)[0]-direct.evaluate(q)[0])<1e-14,"coordinate evaluator differs");
        const auto subset=model.extract({0});
        require(subset.fitDiagnostics().directFitCount==0,"extraction invented a fit for manual coefficients");
        require(std::abs(subset.discreteEigenvalues()[0]-mu)<1e-14,"extraction lost sampling interval");
        path=std::filesystem::temp_directory_path()/("dynsyslab-manual-polynomial-"+std::to_string(
            std::chrono::high_resolution_clock::now().time_since_epoch().count())+".npz");
        core::saveKoopmanModel(path,model);
        const auto loaded=core::loadKoopmanModel(path);
        require(std::abs(loaded.discreteEigenvalues()[0]-mu)<1e-14,"archive lost sampled spectrum");
        require(loaded.spectrum().domain()==core::EigenvalueDomain::discreteTime,"archive lost eigenvalue domain");
        require(std::abs(loaded.evaluateCoordinates(q)[0]-model.evaluateCoordinates(q)[0])<1e-14,"archive changed polynomial");
        std::filesystem::remove(path);
        std::cout<<"Manual polynomial eigenfunctions: PASS\n";return 0;
    } catch(const std::exception& e) {
        if(!path.empty()){std::error_code ignored;std::filesystem::remove(path,ignored);}
        std::cerr<<e.what()<<'\n';return 1;
    }
}
