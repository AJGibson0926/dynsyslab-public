#include <dynsyslab/core/CoordinateClosure.hpp>
#include <dynsyslab/core/System.hpp>

#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace dynsyslab::core;

namespace {
void check(bool condition,const char* message) {
    if(!condition)throw std::runtime_error(message);
}
}

int main() {
    try {
        SystemConfiguration configuration;
        configuration.identifier="rpe";
        configuration.parameters["kappa"]=4.0/3.0;
        const auto system=createSystem(configuration);
        const auto echo=system->observableCatalog().at("echo");
        const auto jet=observableJet(echo,system->driftExpressions(),2);
        check(jet.outputDimension()==3,"scalar order-two jet must contain three coordinates");

        ObservableClosureOptions options;
        options.stateMinimum={0.5,-10.0};
        options.stateMaximum={1.5,10.0};
        options.highestDerivativeOrder=2;
        options.sampleCount=16;
        options.pairCandidateCount=1;
        options.maximumRefinementIterations=20;
        options.collisionTolerance=1.0e-10;
        options.pairSeeds={
            0.7500373016335908,5.627534424684314,
            1.1887083801676752,7.575355314915713};
        const auto report=analyzeObservableClosure(*system,echo,options);
        check(report.orders.size()==3,"closure report derivative orders");
        const auto& witness=report.orders[2].closureObstruction;
        check(witness.found,"known RPE order-two collision was not recovered");
        check(witness.leftCoordinates.size()==3&&witness.successorDerivativeLeft.size()==1,
            "closure witness coordinate dimensions");
        check(witness.normalizedResidual<1.0e-10,"closure witness residual");
        check(std::abs(witness.successorDerivativeLeft[0]-witness.successorDerivativeRight[0])>1000.0,
            "third derivative must separate the known order-two collision");
        std::cout<<"observable jet construction, rank scan, and collision refinement passed\n";
        return 0;
    } catch(const std::exception& error) {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
