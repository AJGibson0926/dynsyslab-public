#include <dynsyslab/core/EigenvalueOptimization.hpp>

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition,const char* message){
    if(!condition)throw std::runtime_error(message);
}

}  // namespace

int main(){
    try{
        constexpr std::size_t samples=24,trajectories=4;
        constexpr double dt=0.1,expected=0.82;
        dynsyslab::core::KoopmanTrajectoryEnsemble ensemble;
        ensemble.samples=samples;ensemble.trajectories=trajectories;
        ensemble.stateDimension=1;
        for(std::size_t sample=0;sample<samples;++sample)
            ensemble.time.push_back(dt*static_cast<double>(sample));
        for(std::size_t trajectory=0;trajectory<trajectories;++trajectory){
            const double amplitude=0.5+static_cast<double>(trajectory);
            for(std::size_t sample=0;sample<samples;++sample)
                ensemble.states.push_back(amplitude*std::pow(
                    expected,static_cast<double>(sample)));
        }
        const auto observable=dynsyslab::core::Observable::stateComponent(
            1,0,"signal","Signal","dimensionless");
        dynsyslab::core::EigenvalueOptimizationOptions options;
        options.initialDiscreteGroups={{{0.65,0.0}}};
        options.lowerBound=-1.05;options.upperBound=1.05;
        options.maximumIterations=200;options.tolerance=1.0e-11;
        const auto result=dynsyslab::core::optimizeTrajectoryEigenvalues(
            ensemble,observable,options);
        require(result.discreteEigenvalues.size()==1,
            "optimizer returned the wrong spectrum size");
        require(std::abs(result.discreteEigenvalues[0]-expected)<1.0e-6,
            "optimizer did not recover the trajectory eigenvalue");
        require(result.diagnostics[0].relativeResidual<1.0e-7,
            "optimized variable-projection residual is too large");
        require(std::abs(std::exp(result.continuousEigenvalues[0]*dt)-expected)<1.0e-10,
            "discrete/continuous conversion is inconsistent");
        std::cout<<"dynsyslab_eigenvalue_optimization_tests: PASS\n";
        return 0;
    }catch(const std::exception& error){
        std::cerr<<"dynsyslab_eigenvalue_optimization_tests: FAIL: "
                 <<error.what()<<'\n';return 1;
    }
}
