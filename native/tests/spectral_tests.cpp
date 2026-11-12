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
    for(std::size_t k=0;k<result.samples;++k) {
        std::complex<double> prediction{};
        for(std::size_t j=0;j<result.rank;++j) {
            prediction+=result.reconstructionMatrix[j]*result.coordinates[k*result.rank+j];
            if(k)require(std::abs(result.coordinates[k*result.rank+j]
                -result.discreteEigenvalues[j]*result.coordinates[(k-1)*result.rank+j])<1e-9,
                "oscillatory modal propagation");
        }
        require(std::abs(prediction-std::cos(.2*k))<1e-9,"oscillatory modal reconstruction");
    }
    dynsyslab::core::SystemConfiguration configuration;
    configuration.identifier="duffing"; configuration.preset="double_well_reference";
    const auto system=dynsyslab::core::createSystem(configuration);
    const auto observable=system->observableCatalog().at("total_energy");
    const auto energy=dynsyslab::core::hankelDynamicModeDecomposition(
        states,80,2,0.01,observable,options);
    require(energy.observableName=="total_energy", "named H-DMD observable");
}
void testJoint() {
    using namespace dynsyslab::core;
    std::vector<double> values;
    for(int k=0;k<100;++k) {
        const double a=std::pow(0.97,k),b=std::pow(0.91,k);
        values.insert(values.end(),{a+b,3*a-2*b,5*b});
    }
    HankelDMDOptions options;options.matrixSize=30;options.samplingWindow=12;
    options.rank=2;options.channelScales={2,7,0.5};options.sampleStride=2;
    for(bool end:{false,true}) {
        options.alignToWindowEnd=end;
        const auto fit=hankelDynamicModeDecomposition(values,100,
            std::vector<std::string>{"x","echo","other"},0.01,options);
        require(fit.rank==2&&fit.matrixRows==36,"joint Hankel geometry");
        require(fit.observableNames.size()==3,"joint names");
        for(std::size_t sample=0;sample<fit.samples;++sample)
            for(std::size_t channel=0;channel<3;++channel) {
                std::complex<double> prediction{};
                for(std::size_t mode=0;mode<fit.rank;++mode)
                    prediction+=fit.reconstructionMatrix[channel*fit.rank+mode]
                        *fit.coordinates[sample*fit.rank+mode];
                require(std::abs(prediction-values[(sample+fit.referenceDelay)*2*3+channel])<1e-10,
                    "joint reconstruction scaling or time alignment");
            }
        for(std::size_t sample=0;sample+1<fit.samples;++sample)
            for(std::size_t mode=0;mode<fit.rank;++mode)
                require(std::abs(fit.coordinates[(sample+1)*fit.rank+mode]
                    -fit.discreteEigenvalues[mode]*fit.coordinates[sample*fit.rank+mode])<1e-9,
                    "joint eigenfunction propagation");
    }
    options.snapshotStride=3;options.samplingWindow=4;options.matrixSize=12;
    const auto strided=hankelDynamicModeDecomposition(values,100,
        std::vector<std::string>{"x","echo","other"},0.01,options);
    require(strided.coordinateSampleOffsets.front()==6
        &&strided.coordinateSampleOffsets[1]==12
        &&strided.coordinateSampleOffsets.back()==74,"snapshot stride alignment");
    for(std::size_t k=1;k<strided.samples;++k)for(std::size_t j=0;j<strided.rank;++j) {
        const double dt=.01*(strided.coordinateSampleOffsets[k]-strided.coordinateSampleOffsets[k-1]);
        require(std::abs(strided.coordinates[k*strided.rank+j]
            -std::exp(strided.continuousEigenvalues[j]*dt)*strided.coordinates[(k-1)*strided.rank+j])<1e-9,
            "snapshot stride confused with operator time step");
    }
    bool rejected=false;
    try{(void)hankelDynamicModeDecomposition(values,100,
        std::vector<std::string>{"x","x","other"},0.01,options);}
    catch(const std::invalid_argument&){rejected=true;}
    require(rejected,"duplicate channels accepted");
}

}
int main() {
    try { testSvd(); testDmd(); testHankel(); testJoint(); std::cout << "dynsyslab_spectral_tests: PASS\n"; return 0; }
    catch (const std::exception& error) { std::cerr << "dynsyslab_spectral_tests: FAIL: " << error.what() << '\n'; return 1; }
}
