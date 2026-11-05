#include <dynsyslab/core/SINDy.hpp>
#include <dynsyslab/core/System.hpp>

#include <cmath>
#include <iostream>
#include <stdexcept>

int main() {
    try {
        constexpr std::size_t samples = 101; const double dt = 0.01;
        std::vector<double> states, velocities, time;
        for (std::size_t index = 0; index < samples; ++index) {
            const double t = index * dt, x = std::exp(-t);
            time.push_back(t); states.push_back(x); velocities.push_back(-x);
        }
        dynsyslab::core::KoopmanSINDyOptions options;
        options.polynomialDegree = 2; options.threshold = 1e-8; options.maximumIterations = 2;
        options.normalizeStates = false; options.materializeTrainingValues = true;
        const auto result = dynsyslab::core::fitKoopmanEigenfunctions(
            states, samples, 1, time, {{-1.0, 0.0}}, velocities, options);
        if (result.featureCount != 3 || result.coefficients.size() != 3)
            throw std::runtime_error("polynomial feature shape");
        const double point[]{0.5};
        const auto value = dynsyslab::core::evaluatePolynomialEigenfunction(result, 0, point);
        const auto gradient = dynsyslab::core::gradientPolynomialEigenfunction(result, 0, point);
        if (!std::isfinite(value.real()) || gradient.size() != 1 || result.residualRms[0] > 1e-8)
            throw std::runtime_error("polynomial eigenfunction fit");
        std::vector<double> dynamicalStates, dynamicalVelocities;
        for (std::size_t index=0;index<samples;++index) {
            const double x=states[index];
            dynamicalStates.insert(dynamicalStates.end(),{x,0.25*x});
            dynamicalVelocities.insert(dynamicalVelocities.end(),{-x,-0.25*x});
        }
        dynsyslab::core::SystemConfiguration configuration;
        configuration.identifier="duffing"; configuration.preset="double_well_reference";
        const auto system=dynsyslab::core::createSystem(configuration);
        const auto coordinate=system->observableCatalog().at("position");
        const auto mapped=dynsyslab::core::fitKoopmanEigenfunctionsInCoordinates(
            dynamicalStates,samples,2,time,{{-1.0,0.0}},coordinate,
            dynamicalVelocities,options);
        if (mapped.stateDimension!=1 || mapped.dynamicalStateDimension!=2
            || mapped.coordinateObservableNames!=std::vector<std::string>{"position"}
            || mapped.residualRms[0]>1e-8)
            throw std::runtime_error("observable-coordinate SINDy fit");
        std::cout << "dynsyslab_sindy_tests: PASS\n"; return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_sindy_tests: FAIL: " << error.what() << '\n'; return 1;
    }
}
