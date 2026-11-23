#pragma once
#include <dynsyslab/math/GeneralizedLaplace.hpp>
#include <dynsyslab/math/HarmonicAnalysis.hpp>
#include <functional>
#include <optional>
#include <string>
#include <vector>
namespace dynsyslab::math {
// Stationary/non-stationary split by asymptotic phase. The isochron through x0 picks the point
// of the cycle that x(t) converges to: Winfree, J. Theor. Biol. 16 (1967) 15-42; Guckenheimer,
// J. Math. Biol. 1 (1975) 259-273; as Koopman eigenfunction level sets, Mauroy, Mezic & Moehlis,
// Physica D 261 (2013) 19-30. The split is that of Susuki, Mauroy & Mezic, SIAM J. Appl. Dyn.
// Syst. 20(4) (2021), eq. (23) and sec. 5.2:
//   y_stationary(t)=sum_m g_m e^{i m (theta0+Omega (t-t0))},  y_transient=y-y_stationary.
struct IsochronSplitOptions {
    std::optional<double> asymptoticPhase; // theta0 from a phase map; empty: estimated from the tail.
    double transientFraction{.6};          // Tail used for the independent phase estimate.
    double noiseFloor{1e-9};               // Relative to the largest transient.
    double phaseTolerance{1e-3};           // Disagreement with the tail that raises a warning.
};
struct IsochronSplit {
    double asymptoticPhase{};
    std::string phaseSource;               // "phase_map" or "tail_estimate".
    std::optional<double> tailPhase;       // Independent estimate from the steady tail.
    std::optional<double> phaseMismatch;   // Wrapped difference modulo the reference symmetry.
    ComplexVector stationary,transient;
    std::optional<double> decayRate;
    double fitResidual{};
    std::size_t strobeCount{};
    std::vector<std::string> warnings;
};
// cycle holds the limit-cycle Fourier coefficients g_m of the observable and Omega.
IsochronSplit isochronSplit(const SampledSignal& trajectory,const FourierSeries& cycle,const IsochronSplitOptions& ={});
using PhaseMap=std::function<double(const std::vector<double>&)>;
// Evaluates a caller-supplied asymptotic phase map at the initial state and splits.
IsochronSplit isochronSplitFromMap(const SampledSignal& trajectory,const FourierSeries& cycle,const PhaseMap& map,
    const std::vector<double>& initialState,IsochronSplitOptions options={});
}
