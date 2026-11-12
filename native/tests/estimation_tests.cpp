#include <dynsyslab/core/Estimation.hpp>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
void require(bool condition,const char* message) {
    if (!condition) throw std::runtime_error(message);
}

dynsyslab::core::KoopmanModel twoModeModel() {
    dynsyslab::core::KoopmanSINDyResult fit;
    fit.continuousEigenvalues={{-1.0,0.0},{-2.0,0.0}};
    fit.coefficients={{1.0,0.0},{0.0,0.0},{0.0,0.0},{1.0,0.0}};
    fit.exponents={1,0,0,1}; fit.center={0.0,0.0}; fit.scale={1.0,1.0};
    fit.stateDimension=2; fit.featureCount=2;
    dynsyslab::core::SystemConfiguration system;
    system.identifier="duffing"; system.preset="double_well_reference";
    dynsyslab::core::KoopmanModel model(std::move(fit),std::move(system));
    model.setReconstruction("x",{{{1.0,0.0},{0.0,0.0}},1,2,0.0,0.0,0.0});
    model.setReconstruction("state",{{{1.0,0.0},{0.0,0.0},
                                       {0.0,0.0},{1.0,0.0}},2,2,0.0,0.0,0.0});
    return model;
}
}

int main() {
    try {
        const std::vector<double> A{0.0,1.0,-2.0,-3.0};
        const std::vector<double> C{1.0,0.0};
        const auto observable=dynsyslab::core::analyzeObservability(A,C,2,1);
        require(observable.rank==2,"observable pair should have full rank");
        require(observable.matrix.size()==4,"observability matrix shape is wrong");
        require(observable.gramian.size()==4,"classical observability Gramian is missing");
        require(observable.method==dynsyslab::core::ObservabilityMethod::classicalMatrix,
                "classical observability method tag is wrong");
        dynsyslab::core::ObservabilityOptions discreteOptions;
        discreteOptions.method=dynsyslab::core::ObservabilityMethod::discreteGramian;
        discreteOptions.blockCount=2;
        const auto discrete=dynsyslab::core::analyzeObservability(
            std::vector<double>{0.5,0.0,0.0,0.25},
            std::vector<double>{1.0,0.0,0.0,1.0},2,2,discreteOptions);
        require(std::abs(discrete.gramian[0]-1.25)<1.0e-12
                && std::abs(discrete.gramian[3]-1.0625)<1.0e-12,
                "finite discrete observability Gramian is wrong");
        dynsyslab::core::ObservabilityOptions finiteOptions;
        finiteOptions.method=dynsyslab::core::ObservabilityMethod::continuousFiniteGramian;
        finiteOptions.horizon=1.0;
        const std::vector<double> stable{-1.0,0.0,0.0,-2.0};
        const std::vector<double> identity{1.0,0.0,0.0,1.0};
        const auto finite=dynsyslab::core::analyzeObservability(
            stable,identity,2,2,finiteOptions);
        require(std::abs(finite.gramian[0]-(1.0-std::exp(-2.0))/2.0)<1.0e-10
                && std::abs(finite.gramian[3]-(1.0-std::exp(-4.0))/4.0)<1.0e-10,
                "finite continuous observability Gramian is wrong");
        dynsyslab::core::ObservabilityOptions infiniteOptions;
        infiniteOptions.method=
            dynsyslab::core::ObservabilityMethod::continuousInfiniteGramian;
        const auto infinite=dynsyslab::core::analyzeObservability(
            stable,identity,2,2,infiniteOptions);
        require(std::abs(infinite.gramian[0]-0.5)<1.0e-10
                && std::abs(infinite.gramian[3]-0.25)<1.0e-10,
                "infinite continuous observability Gramian is wrong");
        const auto coupledInfinite=dynsyslab::core::analyzeObservability(
            A,C,2,1,infiniteOptions);
        require(coupledInfinite.rank==2,
                "coupled infinite Gramian should retain both observable directions");
        const auto oscillatoryInfinite=dynsyslab::core::analyzeObservability(
            std::vector<double>{-1.0,-2.0,2.0,-1.0},C,2,1,infiniteOptions);
        require(oscillatoryInfinite.rank==2,
                "complex-pair infinite Gramian should retain both observable directions");
        require(std::string(dynsyslab::core::observabilityMethodName(infinite.method))
                =="continuous-infinite-gramian","observability method name is wrong");
        bool unstableGramianRejected=false;
        try {
            (void)dynsyslab::core::analyzeObservability(
                std::vector<double>{1.0,0.0,0.0,-2.0},identity,2,2,infiniteOptions);
        } catch (const std::invalid_argument&) { unstableGramianRejected=true; }
        require(unstableGramianRejected,
                "infinite Gramian must reject a non-strictly-stable generator");
        const auto hidden=dynsyslab::core::analyzeObservability(
            std::vector<double>{1.0,0.0,0.0,-2.0},std::vector<double>{0.0,1.0},2,1);
        require(hidden.rank==1,"hidden unstable state should be unobservable");
        bool rejected=false;
        try {
            (void)dynsyslab::core::solveContinuousLQE(
                std::vector<double>{1.0,0.0,0.0,-2.0},std::vector<double>{0.0,1.0},
                std::vector<double>{1.0,0.0,0.0,1.0},std::vector<double>{1.0},2,1,true);
        } catch (const std::invalid_argument&) { rejected=true; }
        require(rejected,"LQE must reject a nondetectable unstable direction");
        const auto lqe=dynsyslab::core::solveContinuousLQE(
            A,C,std::vector<double>{1.0,0.0,0.0,1.0},std::vector<double>{0.5},2,1);
        require(lqe.detectable,"observable pair should be detectable");
        require(lqe.gain.size()==2,"LQE gain shape is wrong");
        require(lqe.relativeResidual<1.0e-6,"LQE CARE residual is too large");
        for (const auto pole:lqe.observerPoles)
            require(pole.real()<0.0,"LQE observer must be stable");
        auto model=twoModeModel();
        const auto design=dynsyslab::core::designKoopmanLQE(
            model,"x","state",std::vector<double>{1.0,0.0,0.0,1.0},
            std::vector<double>{0.25});
        const std::vector<double> time{0.0,0.1,0.2,0.3};
        const std::vector<double> measurements{1.0,0.9,0.8,0.7};
        const auto estimated=dynsyslab::core::runKoopmanLQE(
            model,design,time,measurements,std::vector<double>{0.0,0.0});
        require(estimated.liftedEstimate.size()==8,"KLQE lifted history shape is wrong");
        require(estimated.stateEstimate.size()==8,"KLQE state history shape is wrong");
        require(estimated.measurementEstimate.size()==4,"KLQE measurement history shape is wrong");
        require(estimated.stateEstimate[6]>0.0,"KLQE did not respond to held measurements");
        std::cout<<"native observability and LQE tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr<<"estimation test failure: "<<error.what()<<'\n';
        return 1;
    }
}
