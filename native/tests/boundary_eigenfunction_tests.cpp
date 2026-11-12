#include <dynsyslab/core/BoundaryEigenfunctions.hpp>
#include <dynsyslab/core/KoopmanModel.hpp>

#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>

namespace {
namespace core=dynsyslab::core;
void require(bool condition,const char* message){
    if(!condition)throw std::runtime_error(message);
}
}

int main(){
    try{
        constexpr double dt=0.1;constexpr std::size_t samples=5,trajectories=2;
        core::KoopmanTrajectoryEnsemble ensemble;
        ensemble.samples=samples;ensemble.trajectories=trajectories;
        ensemble.stateDimension=2;ensemble.time={0.0,0.1,0.2,0.3,0.4};
        for(std::size_t trajectory=0;trajectory<trajectories;++trajectory)
            for(std::size_t sample=0;sample<samples;++sample){
                const double slow=std::pow(0.5,static_cast<double>(sample));
                const double x=(trajectory?2.0:1.0)+(trajectory?0.5:2.0)*slow;
                ensemble.states.push_back(x);
                ensemble.states.push_back(static_cast<double>(trajectory)+0.1*sample);
            }
        const auto observable=core::Observable::stateComponent(
            2,0,"x","Position","dimensionless");
        const std::vector<std::complex<double>> eigenvalues{
            {0.0,0.0},{std::log(0.5)/dt,0.0}};
        const auto fit=core::fitBoundaryEigenfunctions(
            eigenvalues,dt,ensemble,observable);
        require(fit.eigenfunctions.samples==samples*trajectories
            &&fit.eigenfunctions.trainingValues.size()==samples*trajectories*2,
            "boundary eigenfunction sample layout is wrong");
        require(fit.relativeResidualByOutput.size()==1
            &&fit.relativeResidualByOutput[0]<1.0e-11,
            "boundary regression did not recover exact modal dynamics");
        core::SystemConfiguration system;
        system.identifier="duffing";system.preset="double_well_reference";
        core::KoopmanModel model(fit.eigenfunctions,std::move(system));
        model.setReconstruction("x",fit.observableReconstruction);
        const auto reconstructed=model.reconstruct("x",ensemble.states.data());
        require(std::abs(reconstructed[0].real()-ensemble.states[0])<1.0e-11,
            "boundary reconstruction C map does not reproduce its observable");
        require(model.eigenfunctionDescriptors()[0].representation
            =="scattered-interpolation",
            "boundary eigenfunction representation was not retained");
        std::cout<<"dynsyslab_boundary_eigenfunction_tests: PASS\n";return 0;
    }catch(const std::exception& error){
        std::cerr<<"dynsyslab_boundary_eigenfunction_tests: FAIL: "
                 <<error.what()<<'\n';return 1;
    }
}
