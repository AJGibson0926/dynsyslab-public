#pragma once
#include <dynsyslab/math/TorusFourier.hpp>
namespace dynsyslab::math {
// Positive finite clock rho(theta)=dt/ds for dtheta/ds=omega.
// h solves omega.grad(h)=rho-mean(rho); eta=theta+omega*h/mean(rho).
// On the supplied polynomial, deta/dt=omega/mean(rho). This also applies to
// any finite circle-cover factor of a solenoid. Unknown tails are excluded.
struct TorusClockTransform {
    TorusFourierSeries clockRate, primitive;
    double meanRate{}, lowerBound{}, upperBound{}, primitiveL1Bound{};
    std::vector<double> straightenedFrequencies;
};
// Requires exact Hermitian coefficient pairs and a strictly positive Fourier
// l1 lower bound (sufficient, not necessary). Nonconstant terms whose divisors
// fall within tolerance are rejected, never discarded. Bounds exclude roundoff.
TorusClockTransform analyzeTorusClock(const TorusFourierSeries&, double divisorTolerance=1e-12);
struct TorusClockSample {
    double time{}, baseElapsed{}, timeResidual{}, clockRate{};
    std::vector<double> phases, straightenedPhases;
    std::vector<std::complex<double>> bareCharacters, correctedCharacters;
};
struct TorusClockEvolution {
    TorusClockTransform transform;
    std::vector<TorusClockSample> samples;
};
// Inverts t=mean(rho)*s+h(theta0+omega*s)-h(theta0) by bounded bisection.
// Times are elapsed physical times, including negative/unordered times.
// The character list uses integer indices on the supplied common torus factor.
TorusClockEvolution evolveTorusClock(const TorusFourierSeries&, const std::vector<double>& initialPhases,
    const std::vector<double>& times, const std::vector<FourierIndex>& characters={},
    double divisorTolerance=1e-12, double relativeTimeTolerance=1e-12);
}
