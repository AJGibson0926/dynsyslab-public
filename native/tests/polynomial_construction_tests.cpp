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
        auto polynomial=core::KoopmanSINDyResult::fromCoefficients(2,powers,coefficients,{lambda});
        const double q[]{3,4};
        require(std::abs(core::evaluatePolynomialEigenfunction(polynomial,0,q)-std::complex<double>(5,4))<1e-14,
            "constant/complex coefficients not preserved");
        require(polynomial.samples==0&&polynomial.directFits==0&&polynomial.residualRms.empty(),"invented fit evidence");
        const auto gradient=core::gradientPolynomialEigenfunction(polynomial,0,q);
        require(gradient==std::vector<std::complex<double>>({{1,0},{0,1}}),"manual analytic gradient failed");
        rejects([&]{(void)polynomial.discreteEigenvalues();});
        polynomial.sampleTime=.1;
        const auto mu=polynomial.discreteEigenvalues()[0];
        require(std::abs(mu-std::exp(lambda*.1))<1e-14,"continuous to discrete failed");
        require(std::abs(polynomial.discreteEigenvalues(.2)[0]-mu*mu)<1e-14,"new dt conversion failed");
        auto sampled=core::KoopmanSINDyResult::fromDiscreteCoefficients(2,powers,coefficients,{mu},.1);
        require(std::abs(sampled.continuousEigenvalues[0]-lambda)<1e-14,"discrete to continuous failed");
        auto aliased=core::KoopmanSINDyResult::fromDiscreteCoefficients(2,powers,coefficients,{mu},.1,{}, {},{1});
        require(std::abs(aliased.discreteEigenvalues()[0]-mu)<1e-13,"branch did not preserve multiplier");
        require(std::abs(aliased.continuousEigenvalues[0].imag()-3-20*std::acos(-1.0))<1e-12,"branch rate incorrect");
        rejects([&]{(void)core::KoopmanSINDyResult::fromDiscreteCoefficients(2,powers,coefficients,{0},.1);});
        rejects([&]{(void)core::KoopmanSINDyResult::fromDiscreteCoefficients(2,powers,coefficients,{mu},0);});
        rejects([&]{(void)core::KoopmanSINDyResult::fromCoefficients(2,powers,{1},{lambda});});
        rejects([&]{(void)core::KoopmanSINDyResult::fromCoefficients(2,powers,{0,0,0},{lambda});});
        rejects([&]{(void)core::KoopmanSINDyResult::fromCoefficients(2,{0,0,0,0},{1,2},{lambda});});
        rejects([&]{(void)core::KoopmanSINDyResult::fromCoefficients(2,powers,coefficients,{lambda},{0,0},{1,0});});
        rejects([&]{(void)core::KoopmanSINDyResult::fromCoefficients(2,powers,coefficients,
            {{std::numeric_limits<double>::quiet_NaN(),0}});});
        core::SystemConfiguration system; system.identifier="duffing";
        const auto direct=core::KoopmanModel::fromPolynomialCoefficients(system,2,powers,coefficients,{lambda});
        require(std::abs(direct.evaluate(q)[0]-std::complex<double>(5,4))<1e-14,"manual KROM evaluation differs");
        sampled.coordinateProvider="causal-history";
        sampled.coordinateObservableNames={"xi_1","xi_2"};
        core::KoopmanModel model(sampled,system);
        require(model.requiresExternalCoordinates(),"history coordinates became physical state");
        require(std::abs(model.evaluateCoordinates(q)[0]-direct.evaluate(q)[0])<1e-14,"coordinate evaluator differs");
        const auto subset=model.extract({0});
        require(subset.eigenfunctions().directFits==0,"extraction invented a fit for manual coefficients");
        require(std::abs(subset.discreteEigenvalues()[0]-mu)<1e-14,"extraction lost sampling interval");
        path=std::filesystem::temp_directory_path()/("dynsyslab-manual-polynomial-"+std::to_string(
            std::chrono::high_resolution_clock::now().time_since_epoch().count())+".npz");
        core::savePolynomialKoopmanModel(path,model);
        const auto loaded=core::loadPolynomialKoopmanModel(path);
        require(std::abs(loaded.discreteEigenvalues()[0]-mu)<1e-14,"archive lost sampled spectrum");
        require(loaded.eigenfunctions().eigenvalueOrigin=="discrete-time","archive lost eigenvalue domain");
        require(std::abs(loaded.evaluateCoordinates(q)[0]-model.evaluateCoordinates(q)[0])<1e-14,"archive changed polynomial");
        std::filesystem::remove(path);
        std::cout<<"Manual polynomial eigenfunctions: PASS\n";return 0;
    } catch(const std::exception& e) {
        if(!path.empty()){std::error_code ignored;std::filesystem::remove(path,ignored);}
        std::cerr<<e.what()<<'\n';return 1;
    }
}
