#include <dynsyslab/core/Extension.hpp>

#include <cmath>
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
        std::cout<<"dynsyslab_extension_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr<<"dynsyslab_extension_tests: FAIL: "<<error.what()<<'\n';
        return 1;
    }
}
