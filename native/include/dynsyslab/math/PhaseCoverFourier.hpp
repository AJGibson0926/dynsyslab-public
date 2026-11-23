#pragma once
#include <dynsyslab/math/TorusFourier.hpp>
#include <cstdint>
#include <string>
namespace dynsyslab::math {
struct PhaseCoverForcingMode {
    FourierIndex index;
    std::string frequencyRatio;
    std::size_t minimumDepth{};
    double detuning{},forcingMagnitude{},correctionMagnitude{};
    bool retained{};
};
struct PhaseCoverForcingResult {
    std::string coverDegree;
    TorusFourierSeries forcing;
    CohomologicalSolution response;
    std::vector<PhaseCoverForcingMode> modes;
};
// Finite characters of Sigma_p x T^d at a supplied common depth D:
// index[0]=k denotes chi_{k/p^D}; subsequent indices are forcing harmonics.
// Composes exact cover labels with the existing finite cohomological solver.
// Retention uses floating detunings and a user tolerance; no KAM certificate,
// exact-resonance proof, or infinite Fourier-tail estimate is returned.
PhaseCoverForcingResult phaseCoverForcing(std::uint32_t prime,std::size_t depth,
    double baseAngularFrequency,const std::vector<double>& forcingFrequencies,
    const std::map<FourierIndex,std::complex<double>>& coefficients,double divisorTolerance);
}
