#pragma once
#include <cstddef>
#include <optional>
#include <vector>
namespace dynsyslab::math {
struct HarmonicResonance {
    std::size_t harmonicA{},harmonicB{};
    double detuning{},relativePhaseDrift{},ratioError{};
    std::optional<double> beatTime;
};
// Tests m*omegaA versus n*omegaB. Numerical equality is not an exact
// rationality certificate. Angles/frequencies use radians and one native clock.
HarmonicResonance harmonicResonance(double omegaA,double omegaB,std::size_t m,
                                  std::size_t n,double observationTime=0);
// Primitive record-small detunings for m,n<=maxOrder, ordered by n.
std::vector<HarmonicResonance> resonanceCandidates(double omegaA,double omegaB,
    std::size_t maxOrder,double observationTime=0);
}
