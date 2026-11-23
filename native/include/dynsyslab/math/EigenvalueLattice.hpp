#pragma once
#include <dynsyslab/math/Resonance.hpp>
#include <complex>
#include <cstddef>
#include <optional>
#include <string>
#include <vector>
namespace dynsyslab::math {
// Near-lattice diagnostics for fitted Koopman spectra. Koopman eigenvalues of
// an autonomous system form a semigroup under products of eigenfunctions:
// Mezic, Nonlinear Dynamics 41 (2005), doi:10.1007/s11071-005-2824-x;
// Budisic, Mohr & Mezic, Chaos 22 (2012), doi:10.1063/1.4772195.
// A match is a numerical CANDIDATE within a declared tolerance. It never
// certifies rationality, an exact eigenfunction product, or a true generator.
enum class EigenvalueClock { Continuous, Discrete };
enum class LatticeRole { Constant, Principal, Combination, Unmatched, Excluded };
struct EigenvalueLatticeOptions {
    EigenvalueClock clock{EigenvalueClock::Continuous};
    // Continuous only: identifies frequencies modulo 2*pi/sampleInterval.
    std::optional<double> sampleInterval;
    std::size_t maximumOrder{4},maximumPrincipals{4};
    // Accepted |detuning| <= max(tolerance, relativeTolerance*|coordinate|),
    // in the working coordinate: rate for continuous, log(mu) for discrete.
    double tolerance{1e-3},relativeTolerance{0},observationTime{0};
    // Explicit generators. Empty selects them greedily by ascending |coordinate|.
    std::vector<std::size_t> principalIndices;
};
struct EigenvalueLatticeEntry {
    std::size_t index{};
    std::complex<double> eigenvalue,coordinate,detuning;
    LatticeRole role{LatticeRole::Unmatched};
    std::vector<unsigned> multiIndex; // Aligned with report.principal.
    unsigned order{};
    double detuningMagnitude{},acceptedTolerance{},phaseDrift{};
    // coordinate = sum m_j s_j + i*aliasTurns*aliasPeriod + detuning.
    long long aliasTurns{};
};
struct EigenvalueLatticeReport {
    std::vector<EigenvalueLatticeEntry> entries;
    std::vector<std::size_t> principal; // Input indices of the generators.
    std::optional<double> aliasPeriod;  // Imaginary-axis period, if any.
    double combinationCount{};
    // Phase-only fraction of the alias circle covered by all tolerance balls.
    // Near one, matches are expected by chance and carry little evidence.
    std::optional<double> chanceCoverage;
    bool rationalityCertified{false};
};
EigenvalueLatticeReport eigenvalueLattice(const std::vector<std::complex<double>>& eigenvalues,
                                          const EigenvalueLatticeOptions& options={});
const char* latticeRoleName(LatticeRole role);

struct FrequencyResonanceOptions {
    std::string source{"spectral_krom"}; // "hankel" or "spectral_krom"; a label only.
    bool radiansPerSample{false};        // Period 2*pi; exclusive with sampleInterval.
    std::optional<double> sampleInterval;
    std::size_t maximumOrder{4},maximumPrincipals{4};
    double tolerance{1e-3},relativeTolerance{0},observationTime{0};
    std::vector<std::size_t> principalIndices;
};
struct PrincipalPairResonance { std::size_t principalA{},principalB{}; HarmonicResonance resonance; };
struct FrequencyResonanceReport {
    std::string source;
    EigenvalueLatticeReport lattice;         // Frequencies folded to i*|omega|.
    std::vector<std::size_t> foldedNegative; // Inputs whose sign was folded.
    std::vector<PrincipalPairResonance> pairs; // Within tolerance, via resonanceCandidates.
};
FrequencyResonanceReport frequencyResonanceReport(const std::vector<double>& frequencies,
                                                  const FrequencyResonanceOptions& options={});
}
