#pragma once
#include <dynsyslab/math/HarmonicAnalysis.hpp>
#include <cstddef>
namespace dynsyslab::math {
// Colbrook, Drysdale & Horning, Rigged DMD, arXiv:2405.00782, Appendix B.
// Equispaced upper-half-plane poles; periodic kernel integrates to one.
// Orders >1 may have negative sidelobes. This is a filter, not a probability.
class RationalCircleKernel {
  public:
    explicit RationalCircleKernel(double bandwidth, std::size_t order=1);
    double operator()(double angularDifference) const;
    double timeMultiplier(long long lag) const;
    // L1 bound for Fourier kernel coefficients at |lag|>retainedLag,
    // including the 1/(2*pi) normalization. Multiplied by ||g|| it bounds
    // unitary spectral-filter truncation error. Exact-arithmetic majorant,
    // not a floating-point certificate or a fitted-model error estimate.
    double symmetricTimeTailBound(std::size_t retainedLag) const;
    double bandwidth() const { return bandwidth_; }
    const std::vector<std::complex<double>>& poles() const { return poles_; }
    const std::vector<std::complex<double>>& residues() const { return residues_; }
  private:
    double bandwidth_;
    std::vector<std::complex<double>> poles_, residues_;
    double decay_,oneMinusDecay_,oneMinusDecaySquared_;
    std::vector<std::complex<double>> poleRotations_;
};
// Each Fourier series is in physical-time phase, not an auxiliary elliptic clock.
FourierSeries spectralPacket(const FourierSeries&, double sampleInterval,
                            double spectralAngle, const RationalCircleKernel&);
struct WeightedOrbitSpectrum { FourierSeries series; double weight{}; };
struct SmoothedOrbitSpectrum {
    std::vector<double> density;
    double totalMass{}, zeroHarmonicMass{};
};
// Weights are explicit quadrature masses, not automatically normalized.
// Finite orbit quadrature does NOT certify whether the limiting spectrum is continuous.
SmoothedOrbitSpectrum smoothOrbitSpectra(const std::vector<WeightedOrbitSpectrum>&,
    double sampleInterval, const std::vector<double>& angles, const RationalCircleKernel&);
struct FrequencyBandPrediction {
    std::vector<std::complex<double>> values;
    std::vector<double> knownCoefficientErrorBound;
};
// Replace each physical angular frequency n*Omega by its nearest multiple of
// bandWidth. Zero width gives exact finite Fourier evolution. Bounds cover only
// supplied coefficients, not unknown Fourier tails or uncertain initial phase.
FrequencyBandPrediction predictFrequencyBands(const FourierSeries&, double initialPhase,
    const std::vector<double>& elapsedTimes, double bandWidth);
}
