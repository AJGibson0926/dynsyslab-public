#include <dynsyslab/core/Spectral.hpp>
#include <dynsyslab/core/System.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
void require(bool value, const char* message) { if (!value) throw std::runtime_error(message); }
void testSvd() {
    const auto result = dynsyslab::core::singularValueDecomposition({3,0,0,2}, 2, 2);
    require(result.rank == 2, "SVD rank");
    require(std::abs(result.singularValues[0] - 3.0) < 1e-12, "SVD leading value");
}
void testDmd() {
    std::vector<double> states;
    for (int sample = 0; sample < 20; ++sample) {
        states.push_back(std::pow(0.9, sample)); states.push_back(std::pow(0.7, sample));
    }
    const auto result = dynsyslab::core::dynamicModeDecomposition(states, 20, 2, 0.1);
    std::vector<double> spectrum;
    for (const auto value : result.discreteEigenvalues) spectrum.push_back(value.real());
    std::sort(spectrum.begin(), spectrum.end());
    require(spectrum.size() == 2, "DMD spectrum size");
    require(std::abs(spectrum[0] - 0.7) < 1e-10 && std::abs(spectrum[1] - 0.9) < 1e-10,
            "DMD discrete eigenvalues");
    require(result.continuousEigenvalues.size() == result.discreteEigenvalues.size(),
            "continuous spectrum conversion");
}
void testHankel() {
    std::vector<double> states;
    for (int sample = 0; sample < 80; ++sample) {
        states.push_back(std::cos(0.2 * sample)); states.push_back(std::sin(0.2 * sample));
    }
    dynsyslab::core::HankelDMDOptions options;
    options.matrixSize = 12; options.samplingWindow = 50; options.observableIndex = 0;
    const auto result = dynsyslab::core::hankelDynamicModeDecomposition(states, 80, 2, 0.01, options);
    require(result.rank >= 2, "Hankel-DMD rank");
    require(result.modalAmplitudes.size() == result.rank, "Hankel modal amplitudes");
    require(result.matrixRows == 50 && result.matrixColumns == 12, "Hankel geometry");
    require(result.coordinates.size()==13*result.rank,"Hankel eigenfunction values");
    require(result.reconstructionMatrix.size()==result.rank,"Hankel reconstruction C");
    dynsyslab::core::SystemConfiguration configuration;
    configuration.identifier="duffing"; configuration.preset="double_well_reference";
    const auto system=dynsyslab::core::createSystem(configuration);
    const auto observable=system->observableCatalog().at("total_energy");
    const auto energy=dynsyslab::core::hankelDynamicModeDecomposition(
        states,80,2,0.01,observable,options);
    require(energy.observableName=="total_energy", "named H-DMD observable");
}
}
int main() {
    try { testSvd(); testDmd(); testHankel(); std::cout << "dynsyslab_spectral_tests: PASS\n"; return 0; }
    catch (const std::exception& error) { std::cerr << "dynsyslab_spectral_tests: FAIL: " << error.what() << '\n'; return 1; }
}
