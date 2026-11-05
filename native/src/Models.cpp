#include "Models.hpp"

#include <cmath>
#include <stdexcept>

namespace dynsyslab::core::detail {
namespace {

constexpr double pi = 3.141592653589793238462643383279502884;

double forcing(double time, const double* parameters, std::size_t offset) {
    return parameters[offset] *
           std::cos(2.0 * pi * parameters[offset + 1] * time + parameters[offset + 2]);
}

}  // namespace

void evaluateRHS(const std::string& model, const double* parameters, double time,
                 const double* state, double* output) {
    const double q = state[0];
    const double v = state[1];
    output[0] = v;
    if (model == "duffing") {
        output[1] = -parameters[0] * v + parameters[1] * q -
                    parameters[2] * q * q * q + forcing(time, parameters, 3);
        return;
    }
    if (model == "van_der_pol") {
        output[1] = parameters[0] * (1.0 - q * q) * v - q + forcing(time, parameters, 1);
        return;
    }
    if (model == "rpe_dimensionless") {
        const double drivePa = forcing(time * parameters[7], parameters, 8);
        const double numerator =
            parameters[1] * std::pow(q, -3.0 * parameters[0]) + parameters[2] -
            parameters[3] - parameters[4] / q - parameters[5] * v / q -
            1.5 * v * v + drivePa / parameters[6];
        output[1] = numerator / q;
        return;
    }
    throw std::invalid_argument("unsupported native ODE model: " + model);
}

std::size_t parameterCount(const std::string& model) {
    if (model == "duffing") return 6;
    if (model == "van_der_pol") return 4;
    if (model == "rpe_dimensionless") return 11;
    throw std::invalid_argument("unsupported native ODE model: " + model);
}

std::size_t stateDimension(const std::string& model) {
    static_cast<void>(parameterCount(model));
    return 2;
}

void validateState(const std::string& model, const double* state) {
    if (!std::isfinite(state[0]) || !std::isfinite(state[1]))
        throw std::runtime_error("native ODE state became non-finite");
    if (std::hypot(state[0], state[1]) > 1.0e8)
        throw std::runtime_error("native ODE state exceeded the built-in validity bound");
    if (model == "rpe_dimensionless" && state[0] <= 0.01)
        throw std::runtime_error("native dimensionless RPE radius left its validity domain");
}

}  // namespace dynsyslab::core::detail
