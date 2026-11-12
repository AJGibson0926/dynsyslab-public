#include <dynsyslab/c_api.h>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

std::string output(dsl_status (*function)(char*, size_t, size_t*)) {
    size_t required = 0;
    require(function(nullptr, 0, &required) == DSL_STATUS_BUFFER_TOO_SMALL,
            "size query status");
    std::vector<char> buffer(required);
    require(function(buffer.data(), buffer.size(), &required) == DSL_STATUS_OK,
            "buffer output status");
    return buffer.data();
}
}

int main() {
    try {
        require(dsl_abi_version() == DSL_ABI_VERSION, "ABI version");
        const auto catalog = output(dsl_system_catalog_json);
        require(catalog.find("duffing") != std::string::npos,
                "ODE catalog JSON");
        require(catalog.find("mackey_glass") != std::string::npos
                && catalog.find("potts") != std::string::npos
                && catalog.find("heat_equation_1d") != std::string::npos
                && catalog.find("\"family\":\"lattice\"") != std::string::npos
                && catalog.find("\"family\":\"pde\"") != std::string::npos,
                "non-ODE catalog JSON");
        dsl_system_handle* system = nullptr;
        require(dsl_system_create_json(
            "{\"identifier\":\"rpe\",\"preset\":\"air_water_20c\"}", &system)
            == DSL_STATUS_OK, "create RPE");
        size_t dimension = 0;
        require(dsl_system_dimension(system, &dimension) == DSL_STATUS_OK && dimension == 2,
                "RPE dimension");
        const double state[]{1.0, 0.0}; double derivative[2]{};
        require(dsl_system_rhs(system, 0.0, state, 2, 0.0, derivative, 2)
            == DSL_STATUS_OK, "RPE RHS");
        require(std::abs(derivative[1]) < 1.0e-10, "RPE equilibrium");
        size_t required = 0;
        require(dsl_system_configuration_json(system, nullptr, 0, &required)
            == DSL_STATUS_BUFFER_TOO_SMALL, "resolved configuration size");
        std::vector<char> configuration(required);
        require(dsl_system_configuration_json(system, configuration.data(),
            configuration.size(), &required) == DSL_STATUS_OK, "resolved configuration");
        const std::string text(configuration.data());
        require(text.find("equations_fingerprint") != std::string::npos,
                "equations fingerprint");
        require(text.find("characteristic_time_s") != std::string::npos,
                "dimensionless scaling");
        dsl_system_destroy(system);
        std::cout << "dynsyslab_c_api_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_c_api_tests: FAIL: " << error.what() << '\n';
        return 1;
    }
}
