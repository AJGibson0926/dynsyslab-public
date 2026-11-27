#include <dynsyslab/core/Extension.hpp>

#include <cmath>
#include <cstdint>
#include <vector>
#include <iostream>
#include <stdexcept>

namespace {
void require(bool condition,const char* message) {
    if (!condition) throw std::runtime_error(message);
}
}

int main(int count,char** values) {
    try {
        require(count==2,"extension library path argument");
        const auto builtinCount=dynsyslab::core::availableSystems().size();
        auto& registry=dynsyslab::core::ExtensionRegistry::instance();
        require(registry.load(values[1])=="dynsyslab.test","extension identifier");
        require(registry.load(values[1])=="dynsyslab.test","idempotent extension load");
        require(registry.extensions().size()==1,"extension inventory");
        require(dynsyslab::core::availableSystems().size()==builtinCount+1,
                "extension system catalog composition");
        dynsyslab::core::SystemConfiguration configuration;
        configuration.identifier="test_oscillator";
        auto system=dynsyslab::core::createSystem(configuration);
        require(system->descriptor().stateDimension==2,"extension system dimension");
        const double state[]{1.0,0.0}; double derivative[2]{};
        system->rhsWithInput(0.0,state,0.0,derivative);
        require(std::abs(derivative[1]+4.0)<1.0e-14,"extension RHS");
        require(system->configuration().parameters.at("omega")==2.0,
                "extension resolved parameter");
        require(system->observableCatalog().at("energy").evaluateReal(state)[0]==2.0,
                "extension observable");
        require(registry.execute("dynsyslab.test","test","{}").find("true")
            !=std::string::npos,"extension generic operation");
        // Delay, field and lattice entry points.
        require(registry.models().size()==3,"extension model inventory");
        const auto delay=registry.model("test_delay");
        require(delay&&delay->family=="dde","extension delay model");
        const double window[]{0.5,2.0}; double rate{};
        delay->delayRhs(window,2,&rate,1);
        require(rate==-2.0,"extension delay rhs");
        const auto field=registry.model("test_field");
        require(field&&field->family=="pde","extension field model");
        const auto u=field->initialField("{}",8);
        std::vector<double> du(8);
        field->fieldRhs(0.0,u.data(),8,du.data());
        require(std::abs(du[2]-64.0*(2.0*std::sin(0.7853981633974483)-2.0))<1.0e-9,"extension field rhs");
        const auto lattice=registry.model("test_lattice");
        require(lattice&&lattice->family=="lattice","extension lattice model");
        std::int32_t sites[4]{}; double mean{};
        lattice->latticeInitial(1,sites,4);
        lattice->latticeSweep(1,0,sites,4,&mean,1);
        require(mean==-1.0&&sites[0]==-1,"extension lattice sweep");
        std::cout<<"dynsyslab_extension_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr<<"dynsyslab_extension_tests: FAIL: "<<error.what()<<'\n';
        return 1;
    }
}
