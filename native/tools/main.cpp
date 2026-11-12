#include <dynsyslab/core/ODE.hpp>
#include <dynsyslab/core/Version.hpp>

#include <exception>
#include <iomanip>
#include <iostream>
#include <string>

namespace {

void capabilities() {
    std::cout << "{\"schema\":\"dynsyslab.native-capabilities/v1\",\"version\":\""
              << dynsyslab::core::version
              << "\",\"python_required\":false,\"fixed_step_models\":[";
    const auto models = dynsyslab::core::fixedStepModels();
    for (std::size_t index = 0; index < models.size(); ++index) {
        if (index) std::cout << ',';
        std::cout << '\"' << models[index] << '\"';
    }
    std::cout << "]}\n";
}

int solveDuffing() {
    dynsyslab::core::FixedIntegrationRequest request;
    request.model = "duffing";
    request.parameters = {0.2, 1.0, 1.0, 0.0, 0.0, 0.0};
    request.initialStates = {0.1, 0.0};
    request.finalTime = 0.01;
    request.stepSize = 0.001;
    const auto result = dynsyslab::core::integrateFixed(request);
    const auto* final = result.state(0, result.samples - 1);
    std::cout << std::setprecision(17)
              << "{\"schema\":\"dynsyslab.native-solve/v1\",\"samples\":" << result.samples
              << ",\"final_state\":[" << final[0] << ',' << final[1] << "]}\n";
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const std::string command = argc > 1 ? argv[1] : "capabilities";
        if (command == "capabilities") capabilities();
        else if (command == "self-test") return solveDuffing();
        else {
            std::cerr << "usage: dynsyslab_native_cli [capabilities|self-test]\n";
            return 2;
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
