#pragma once
#include <dynsyslab/math/LaplaceTransform.hpp>
#include <map>
#include <vector>
namespace dynsyslab::math {
// Generalized Laplace Analysis: Mezic, Nonlinear Dyn. 41 (2005), doi:10.1007/s11071-005-2824-x;
// Budisic, Mohr & Mezic, Chaos 22 (2012) 047510, doi:10.1063/1.4772195. Modes are projected from
// slowest decay (largest Re lambda) to fastest, subtracting each before the next, because the
// weight e^{-lambda t} amplifies everything slower than lambda.
//
// Time averages use the C-infinity bump window of Das & Yorke, Nonlinearity 31 (2018) 491-501,
// in place of the plain Cesaro mean. An optional direct solve of the coupled
// window equations (a numerical refinement beyond the cited papers) removes cross-leakage.
struct GeneralizedLaplaceOptions {
    bool refine{true};            // Solve the coupled finite-window equations after the ordered pass.
    double maximumCondition{1e12}; // Keep the ordered-pass values above this conditioning.
    // Each mode's window ends where |e^{-lambda t}| reaches this limit.
    double amplificationLimit{1e6};
};
struct GeneralizedLaplaceMode {
    std::complex<double> eigenvalue,coefficient;  // y(t) ~ sum coefficient * e^{lambda (t-t0)}.
    std::complex<double> sequentialCoefficient;   // Classical ordered GLA, before refinement.
    double horizon{};                             // Averaging window length used.
};
struct GeneralizedLaplaceResult {
    std::vector<GeneralizedLaplaceMode> modes;    // Ordered slow to fast.
    ComplexVector residual;
    double residualRms{},signalRms{};
    double refinementCondition{};
    bool refined{false};
};
GeneralizedLaplaceResult generalizedLaplaceAnalysis(const SampledSignal&,ComplexVector eigenvalues,
    const GeneralizedLaplaceOptions& ={});

// Fourier coefficients c_n in y ~ sum c_n e^{i n w t}, from a bump-weighted average over the
// trailing (1-transientFraction) of the record, trimmed to whole periods.
std::map<int,std::complex<double>> estimateHarmonics(const SampledSignal&,double frequency,int maximumHarmonic,
    double transientFraction=.5);

// Stationary/non-stationary split of a limit-cycle trajectory (Susuki, Mauroy & Mezic 2021,
// eq. (23) and sec. 5.2). stationary is the harmonic series of the tail; transient is the rest.
// decayRate is a log-linear fit of |transient| on stroboscopic samples t0+kT, the check of the
// Floquet exponent that the paper makes in Figure 1.
struct StationaryDecomposition {
    std::map<int,std::complex<double>> harmonics;
    ComplexVector stationary,transient;
    double decayRate{},fitResidual{};
    std::size_t strobeCount{};
    bool decayRateAvailable{false};
};
StationaryDecomposition decomposeStationary(const SampledSignal&,double frequency,int maximumHarmonic,
    double transientFraction=.5,double noiseFloor=1e-10);

// Order probe at lambda from data. A simple pole has a steady weighted mean of e^{-lambda t}y;
// a double pole grows linearly, so the second-half to first-half ratio tends to 3.
struct PoleOrderProbe {double ratio{};std::size_t order{};};
PoleOrderProbe probePoleOrder(const SampledSignal&,std::complex<double> lambda);
}
