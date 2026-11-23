#pragma once
#include <dynsyslab/math/HarmonicAnalysis.hpp>
#include <string>
namespace dynsyslab::math {
struct AnalyticFourierExpansion {
    FourierSeries series;
    // Bound for the omitted infinite series in exact arithmetic; not a
    // floating-point/interval enclosure and not a bound for rounded inputs.
    double uniformTailBound{};
    double parameterPeriod{};
};
// Fourier variable phi=2*pi*u/P; sn/cn P=4K(m), dn P=2K(m).
// The declared physical angularFrequency is independent of this u period.
AnalyticFourierExpansion jacobiFourier(const std::string& kind,double parameter,
    int maximumHarmonic,double angularFrequency=1);
// log(1+k*cos(phi)), |k|<1, with a closed-form coefficient and tail formula.
AnalyticFourierExpansion logarithmicCosineFourier(double k,int maximumHarmonic,
    double angularFrequency=1);
// cos(E) or sin(E) as functions of mean anomaly M=E-e*sin(E), 0<=e<1.
// Classical Fourier-Bessel coefficients, NOT samples at uniform eccentric anomaly.
// No infinite-tail or floating-point bound is supplied by this API.
FourierSeries keplerAnomalyFourier(bool sine,double eccentricity,
    int maximumHarmonic,double angularFrequency=1);
}
