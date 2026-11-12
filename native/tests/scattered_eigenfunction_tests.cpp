#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/KoopmanModel.hpp>

#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace {
namespace core = dynsyslab::core;

void require(bool condition,const char* message) {
    if(!condition)throw std::runtime_error(message);
}

core::KoopmanModel model() {
    core::ScatteredEigenfunctionData functions;
    functions.continuousEigenvalues={{-0.1,1.0},{-0.1,-1.0}};
    functions.trainingStates={-1.0,-1.0,-1.0,1.0,1.0,-1.0,1.0,1.0};
    functions.trainingValues={
        {-3.0,0.0},{1.0,0.0},
        {1.0,0.0},{-3.0,0.0},
        {-1.0,0.0},{3.0,0.0},
        {3.0,0.0},{-1.0,0.0}};
    functions.samples=4;functions.stateDimension=2;functions.neighbors=4;
    functions.provenance="native scattered test";
    core::SystemConfiguration system;
    system.identifier="duffing";system.preset="double_well_reference";
    return core::KoopmanModel(std::move(functions),std::move(system));
}

}  // namespace

int main() {
    std::filesystem::path path;
    try {
        auto source=model();
        const double trainingPoint[]{1.0,1.0};
        const auto exact=source.evaluate(trainingPoint);
        require(exact==std::vector<std::complex<double>>{{3.0,0.0},{-1.0,0.0}},
            "scattered interpolant did not preserve training values");
        const double interior[]{0.25,-0.4};
        const auto input=source.exactControlMap(interior);
        require(input.size()==2&&std::isfinite(input[0].real())
            &&std::isfinite(input[1].real()),
            "finite-difference scattered B_Psi is unavailable");
        require(source.eigenfunctionRepresentation()=="scattered-interpolation"
            &&source.eigenfunctionDescriptors()[0].gradientBackend=="finite-difference",
            "eigenfunction representation descriptor is wrong");

        path=std::filesystem::temp_directory_path()/(
            "dynsyslab-scattered-krom-"+std::to_string(
                std::chrono::high_resolution_clock::now().time_since_epoch().count())+".npz");
        core::saveKoopmanModel(path,source,
            R"({"format":"dynsyslab.derived/v1","label":"Scattered KROM"})");
        const auto loaded=core::loadKoopmanModel(path);
        require(loaded.evaluate(trainingPoint)==exact,
            "scattered KROM archive did not round trip");
        const auto reduced=loaded.extract({1});
        require(reduced.liftedDimension()==1
            &&reduced.evaluate(trainingPoint)[0]==exact[1]
            &&reduced.scatteredEigenfunctions()->trainingValues.size()==4,
            "scattered KROM extraction did not create an independent subset");
        std::error_code ignored;std::filesystem::remove(path,ignored);
        std::cout<<"dynsyslab_scattered_eigenfunction_tests: PASS\n";return 0;
    } catch(const std::exception& error) {
        if(!path.empty()){std::error_code ignored;std::filesystem::remove(path,ignored);}
        std::cerr<<"dynsyslab_scattered_eigenfunction_tests: FAIL: "
                 <<error.what()<<'\n';return 1;
    }
}
