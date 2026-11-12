#include <dynsyslab/core/HistoryKoopman.hpp>

#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

}  // namespace

int main() {
    try {
        constexpr std::size_t samples = 640;
        constexpr std::size_t delays = 80;
        constexpr double step = 0.01;
        constexpr double frequency = 5.0;
        std::vector<std::vector<double>> echoes(3);
        for (std::size_t run = 0; run < echoes.size(); ++run) {
            echoes[run].resize(samples);
            const double phase = 0.6 * static_cast<double>(run);
            for (std::size_t sample = 0; sample < samples; ++sample)
                echoes[run][sample] = (1.0 + 0.2 * static_cast<double>(run))
                    * std::sin(frequency * step * static_cast<double>(sample) + phase);
        }
        dynsyslab::core::PooledHankelOptions hankelOptions;
        hankelOptions.delays = delays;
        hankelOptions.rank = 2;
        hankelOptions.oversampling = 2;
        hankelOptions.subspaceIterations = 6;
        const auto basis = dynsyslab::core::fitPooledHankelBasis(echoes, hankelOptions);
        require(basis.trainingTrajectories == echoes.size(),
                "pooled Hankel trajectory count changed");
        require(basis.trainingWindows == echoes.size() * (samples - delays + 1),
                "pooled Hankel crossed a trajectory boundary");
        require(basis.cumulativeEnergy.back() > 1.0 - 1.0e-9,
                "rank-two sinusoid did not retain its Hankel energy");

        const auto offline = dynsyslab::core::projectHistory(echoes[1], basis);
        dynsyslab::core::CausalHistoryProjector online(basis);
        std::size_t row = 0;
        for (const double observation : echoes[1]) {
            const auto coordinate = online.push(observation);
            if (!coordinate) continue;
            require(coordinate->size() == basis.rank,
                    "online history coordinate rank changed");
            for (std::size_t mode = 0; mode < basis.rank; ++mode)
                require(std::abs((*coordinate)[mode]
                    - offline.coordinates[row * basis.rank + mode]) < 2.0e-10,
                    "online history projection differs from offline projection");
            ++row;
        }
        require(row == offline.windows, "online history projection window count changed");

        std::vector<std::vector<double>> coordinateRuns(4);
        for (std::size_t run = 0; run < coordinateRuns.size(); ++run) {
            constexpr std::size_t count = 900;
            coordinateRuns[run].resize(count * 2);
            const double phase = 0.4 * static_cast<double>(run);
            const double amplitude = 0.8 + 0.15 * static_cast<double>(run);
            for (std::size_t sample = 0; sample < count; ++sample) {
                const double time = step * static_cast<double>(sample);
                const double envelope = amplitude * std::exp(-0.08 * time);
                coordinateRuns[run][sample * 2] = envelope
                    * std::cos(frequency * time + phase);
                coordinateRuns[run][sample * 2 + 1] = envelope
                    * std::sin(frequency * time + phase);
            }
        }
        dynsyslab::core::PooledPolynomialEDMDOptions edmdOptions;
        edmdOptions.degree = 1;
        edmdOptions.pairStride = 2;
        edmdOptions.sampleTime = step;
        const auto edmd = dynsyslab::core::fitPooledPolynomialEDMD(
            coordinateRuns, 2, edmdOptions);
        dynsyslab::core::EDMDModeSearchOptions search;
        search.minimumContinuousImaginary = 4.5;
        search.maximumContinuousImaginary = 5.5;
        const auto candidates = dynsyslab::core::selectPolynomialEDMDModes(
            edmd, coordinateRuns, search);
        require(!candidates.empty(), "native pooled EDMD found no rotating mode");
        require(std::abs(candidates.front().continuousEigenvalue.real() + 0.08) < 1.0e-7,
                "native pooled EDMD decay rate changed");
        require(std::abs(candidates.front().continuousEigenvalue.imag() - frequency) < 1.0e-7,
                "native pooled EDMD frequency changed");
        require(candidates.front().meanPairResidual < 1.0e-9,
                "native pooled EDMD one-pair residual is too large");
        const auto canonical = dynsyslab::core::selectedPolynomialEigenfunction(
            edmd, candidates.front());
        require(canonical.stateDimension == 2 && canonical.featureCount == 3
            && canonical.continuousEigenvalues.size() == 1,
                "selected native eigenfunction did not use the canonical contract");
        std::cout << "native causal-history and pooled Koopman tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "history Koopman test failure: " << error.what() << '\n';
        return 1;
    }
}
