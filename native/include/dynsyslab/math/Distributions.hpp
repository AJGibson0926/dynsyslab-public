#pragma once
#include <complex>
#include <vector>
namespace dynsyslab::math {
// Sum_j coefficients[j]*delta^(j)(x-location). Pair with actual derivatives
// of a smooth test function at location, NOT factorial-scaled Taylor terms.
struct PointDistribution {
    double location{};
    std::vector<std::complex<double>> coefficients;
};
std::complex<double> pairTestJet(const PointDistribution& distribution,
                               const std::vector<std::complex<double>>& testDerivatives);
PointDistribution distributionDerivative(const PointDistribution& distribution,unsigned order=1);
PointDistribution multiplySmooth(const PointDistribution& distribution,
                                 const std::vector<std::complex<double>>& smoothDerivatives);
// Finite quadrature approximation to <e^(i*n*phi) dt/T, test>. Pairing is
// bilinear. With normalized=false, dt is used instead of dt/T.
std::complex<double> periodicOrbitPairing(const std::vector<std::complex<double>>& testSamples,
    int harmonic,double period,bool normalized=true);
// Coarea quadrature for delta(H-E). Positive gradient norms are mandatory.
std::complex<double> levelSetDeltaPairing(const std::vector<std::complex<double>>& testValues,
    const std::vector<double>& surfaceWeights,const std::vector<double>& gradientNorms);
struct SpectralDensityRoot {
    double sourceDensity{},frequencyDerivative{};
    int harmonic{};
    std::complex<double> observableCoefficient;
};
double spectralPushforwardDensity(const std::vector<SpectralDensityRoot>& regularRoots);
// An explicitly chosen smooth approximation, never a pointwise value of delta.
double gaussianDeltaDerivative(double x,double epsilon,unsigned order=0);
}
