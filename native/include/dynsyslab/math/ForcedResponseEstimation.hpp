#pragma once
#include <dynsyslab/math/GeneralizedLaplace.hpp>
#include <map>
#include <optional>
#include <string>
#include <vector>
namespace dynsyslab::math {
// H_n from forced steady-state data. Definition 4 and Theorem 1 of Susuki, Katayama, Mauroy &
// Mezic, arXiv:2603.05771v2: under complex u=u0 e^{iwt}, y_n(t)=H_n (u0 e^{iwt})^n.
struct ForcedHarmonicEstimate {
    std::map<int,std::complex<double>> coefficients;  // c_n of e^{i n w t}.
    std::map<int,std::complex<double>> response;      // H_n=c_n/u0^n, n>=1; complex input only.
    double tailRelativeResidual{};                    // Unexplained part of the steady tail.
    std::vector<std::string> warnings;
};
ForcedHarmonicEstimate estimateForcedHarmonics(const SampledSignal& output,double frequency,double amplitude,
    int maximumHarmonic,bool complexInput=true,double transientFraction=.5);
// EDMD (Williams, Kevrekidis & Rowley, J. Nonlinear Sci. 25 (2015) 1307-1346) on a monomial
// dictionary of the skew-product state (x,u). Since u^n is an exact eigenfunction, H_n is the
// output entry of the Koopman mode at mu=e^{i n w h} normalized to 1 on the u^n coordinate.
struct SkewDmdInput {
    std::vector<ComplexVector> states;   // Snapshot rows, each of the state dimension.
    ComplexVector inputs;                // u at the same uniformly spaced times.
    double sampleInterval{},frequency{};
    std::size_t outputState{0};
    int maximumHarmonic{2},maximumDegree{2};
    std::vector<std::vector<int>> dictionary;  // Exponents of (x...,u); empty uses all up to degree.
    double rankTolerance{1e-11};
};
struct SkewDmdHarmonic {
    int harmonic{};
    std::optional<std::complex<double>> response;
    double modeResidual{};          // |mu-(K v)_{u^n}| for the normalized mode.
    double spectrumDistance{};      // Distance from mu to the nearest eigenvalue of K.
};
struct SkewDmdResult {
    bool feasible{false};
    std::size_t dictionarySize{},rank{},snapshots{};
    double closureResidual{},condition{};   // ||Psi+ - K Psi||_F/||Psi+||_F.
    std::vector<SkewDmdHarmonic> harmonics;
    std::vector<std::string> warnings;
};
SkewDmdResult skewProductDmd(const SkewDmdInput&);
// Subharmonic response on the n-fold cover of the forcing circle. u^{1/n} is not single-valued
// on the state space, so H_{1/n} (Definition 5) is meaningful only up to an n-th root of unity:
// a shift by one forcing period multiplies it by e^{2 pi i/n}. Reported invariants are |H| and
// H^n. Cover labels come from phaseCoverForcing when n is a prime power.
struct SubharmonicLabel {int index{};std::string frequencyRatio;std::size_t minimumDepth{};};
struct SubharmonicEstimate {
    int cover{};
    std::map<int,std::complex<double>> coefficients;  // c_k of e^{i k (w/n) t}.
    double subharmonicFraction{};                     // AC energy at k not divisible by n.
    bool detected{false},coverLabelled{false};
    std::optional<double> magnitude,phaseModulo;      // |H| and arg(H) mod 2 pi/n.
    std::optional<std::complex<double>> invariantPower; // H^n.
    ComplexVector branches;                           // All n admissible values of H_{1/n}.
    std::vector<SubharmonicLabel> labels;
    std::vector<std::string> warnings;
};
SubharmonicEstimate subharmonicResponse(const SampledSignal& output,double frequency,double amplitude,int cover,
    int maximumIndex,double transientFraction=.5,double detectionTolerance=1e-6);
}
