#pragma once
#include <dynsyslab/math/HarmonicAnalysis.hpp>
namespace dynsyslab::math {
// g(alpha), w(alpha)=dt/ds>0, alpha=auxiliaryFrequency*s. Physical phase is
// theta=alpha+beta(alpha), beta(0)=0, beta'=w/mean(w)-1. All inputs are finite
// Fourier polynomials. No elliptic-function or differential-equation assumption.
struct PeriodicClockFourier {
    std::vector<FourierSeries> channels;
    FourierSeries phaseCorrection;
    double meanClockRate{}, physicalPeriod{}, minimumSampledClockRate{};
    // A sufficient lower bound for the supplied clock polynomial in exact
    // arithmetic, from Fourier l1 and sample/Lipschitz bounds; not interval-certified.
    double clockPolynomialLowerBound{};
};
// Periodic trapezoidal quadrature in auxiliary phase, with the clock Jacobian.
// Refine samples to assess quadrature; neither omitted input-series tails nor
// output-series tails are bounded by this routine. Frequencies must agree.
PeriodicClockFourier reparameterizePeriodicFourier(const std::vector<FourierSeries>& channels,
    const FourierSeries& clockRate, std::size_t samples, int maximumHarmonic);
struct ClockCoefficientExpansion {
    std::complex<double> coefficient;
    double exactArithmeticRemainderBound{};
};
// Sum_{k=0}^order (-i*n)^k/k! [g*(w/w0)*beta^k]_n, by finite Fourier
// convolutions. Entire exponential series; the remainder bound excludes roundoff
// and omitted terms in the supplied g and w. Large n*||beta|| causes cancellation.
ClockCoefficientExpansion periodicClockCoefficientExpansion(const FourierSeries& channel,
    const FourierSeries& clockRate, int harmonic, unsigned order);
}
