#pragma once
#include <dynsyslab/core/PeriodicAnalysis.hpp>
namespace dynsyslab::core {
struct CoefficientShiftSample {
    double shift{}, relativeFrequencyError{}, initialStateError{}, closureError{};
    std::vector<double> maximumCoefficientError, coefficientL2Error, referenceL2Norm;
};
struct CoefficientShiftCheck {
    double period{}, angularFrequency{}, closureError{};
    std::size_t channels{}, samples{};
    int maximumHarmonic{};
    std::vector<CoefficientShiftSample> shifts;
};
// Numerical Fourier time-translation check for an AUTONOMOUS periodic flow:
// c_n(Phi^s x)=exp(i*n*omega*s)c_n(x). Each shifted state is independently
// re-instantiated as an IC. Native FFT; no fitted model or phase inverse.
// This checks the supplied family/observable implementation, not a theorem
// or the omitted Fourier tail. Degenerate time frequencies are already summed.
CoefficientShiftCheck checkPeriodicCoefficientShift(
    const AnalyticSolutionFamily&,const AnalyticState& initial,
    const std::vector<Observable>& observables,std::size_t samples,
    int maximumHarmonic,const std::vector<double>& shifts,
    double initialTime=0,double explicitPeriod=0);
}
