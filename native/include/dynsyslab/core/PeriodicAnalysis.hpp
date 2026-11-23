#pragma once
#include <dynsyslab/core/AnalyticSolution.hpp>
#include <dynsyslab/core/Observable.hpp>
namespace dynsyslab::core {
struct PeriodicSamples {
    double period{}, angularFrequency{}, closureError{};
    std::size_t stateDimension{}, channels{};
    std::vector<double> times, states;
    std::vector<std::complex<double>> observables; // sample-major
};
// Uniform physical/native-time samples on [t0,t0+T), no repeated endpoint.
// An explicit period overrides trajectory metadata; endpoint agreement is a
// diagnostic, not a proof of periodicity or of the fundamental period.
PeriodicSamples samplePeriodicTrajectory(const AnalyticTrajectory&, std::size_t stateDimension,
    const std::vector<Observable>&, std::size_t samples, double initialTime=0,
    double explicitPeriod=0);
}
