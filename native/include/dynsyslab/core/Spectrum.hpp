#pragma once
// The portion of a Koopman spectrum a method computed, with the origin (method, settings, inputs) that
// produced it. By Lebesgue decomposition (SpectralDecomposition.hpp) it has
//   an atomic part   the eigenpairs: each eigenvalue with its stable mode id, its residual when the
//                    method computed one (Eigenpair.hpp), and the eigenfunctions as named
//                    representations (one output per pair) with their values at the training samples;
//   a continuous part bands: subsets of C on a grid, carrying no eigenvalue (empty unless a
//                    method estimated one: Rigged DMD, the autocorrelation Fejer estimate);
//   observables      the decomposition of each observable over both parts: its Koopman modes (one
//                    coefficient per pair) and its measure on each spectral band, with the
//                    band's generalized eigenfunctions for it (wave packets);
// and optional evidence and the operator itself:
//   pseudospectrum   the minimum residual on a grid of points z, with its minimizers
//                    (Pseudospectrum.hpp; Trefethen & Embree 2005; Colbrook & Townsend 2024);
//   compression      the compressed operator the pairs were computed from, when the method kept one
//                    (KoopmanCompression.hpp): what Rigged DMD and operator rollout act on.
// Selecting pairs keeps their residuals, ids, eigenfunctions and Koopman modes; the continuous part,
// the measures, the pseudospectrum and the operator describe the operator, not the selection, and are
// kept whole.
//
// Each mode stores its value in its own time domain (a continuous rate lambda, or a sampled
// multiplier mu = exp(lambda dt)) plus the selected logarithm branch, so conversion is explicit and
// lossless. A zero sampled multiplier is valid but has no finite continuous rate. `eigenvalue(j)`
// and `eigenvalues()` hand out Eigenvalue objects, with conjugate partners paired.
#include <dynsyslab/core/Eigenpair.hpp>
#include <dynsyslab/core/Pseudospectrum.hpp>
#include <dynsyslab/core/SpectralDecomposition.hpp>
#include <dynsyslab/core/TimeSeriesRecord.hpp>

#include <complex>
#include <cstddef>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace dynsyslab::core {
struct KoopmanCompression;
class ObservableCatalog;
enum class EigenvalueDomain { continuousTime, discreteTime };

class Spectrum {
  public:
    Spectrum() = default;
    static Spectrum fromContinuousRates(std::vector<std::complex<double>> rates, double sampleInterval = 0);
    static Spectrum fromDiscreteMultipliers(std::vector<std::complex<double>> multipliers, double sampleInterval,
                                            std::vector<int> logarithmBranches = {});
    // Sampled multipliers together with the rates the method chose for them (branches inferred
    // from the rates; both values are kept exactly).
    static Spectrum fromDiscrete(std::vector<std::complex<double>> multipliers,
                                 const std::vector<std::complex<double>>& rates, double sampleInterval,
                                 DataOrigin origin = {});
    // Compatibility import: preserve archived rates exactly, including their selected branches,
    // rather than round-tripping through exp/log.
    static Spectrum fromStoredRates(std::vector<std::complex<double>> rates, double sampleInterval,
                                    EigenvalueDomain origin, std::vector<int> logarithmBranches = {});
    // Lossless archive round-trip: validate and retain both recorded values.
    static Spectrum fromStoredDiscrete(std::vector<std::complex<double>> multipliers, double sampleInterval,
                                       std::vector<int> logarithmBranches,
                                       const std::vector<std::complex<double>>& rates);

    EigenvalueDomain domain() const noexcept { return domain_; }
    double sampleInterval() const noexcept { return sampleInterval_; }
    std::size_t size() const noexcept { return values_.size(); }
    const DataOrigin& origin() const noexcept { return origin_; }
    Spectrum& withOrigin(DataOrigin origin) { origin_ = std::move(origin); return *this; }
    // Values in the stored domain (rates or multipliers).
    const std::vector<std::complex<double>>& values() const noexcept { return values_; }
    const std::vector<int>& logarithmBranches() const noexcept { return branches_; }
    bool hasContinuousRates() const noexcept;
    std::optional<std::complex<double>> continuousRate(std::size_t mode) const;
    const std::vector<std::complex<double>>& continuousRates() const;
    std::vector<std::complex<double>> discreteMultipliers(double sampleInterval = 0) const;
    // Nonzero sampled modes use the recorded logarithm branch between samples.
    // Zero multipliers only admit nonnegative integer sampling intervals.
    std::complex<double> transition(std::size_t mode, double elapsed) const;
    // The pairs at positions `modes`, in that order, with their residuals, stable mode ids,
    // eigenfunctions and Koopman modes; the continuous part, the measures (whose atoms are derived from
    // the Koopman modes), the pseudospectrum and the operator are kept whole.
    Spectrum extract(const std::vector<std::size_t>& modes) const;
    // Stable mode ids. A fit numbers its pairs 0..n-1, and a pair keeps its id through every selection,
    // so references to pairs (a child krom's parent link, atoms) name ids, not positions.
    std::size_t modeId(std::size_t mode) const;
    std::vector<std::size_t> modeIds() const;
    std::optional<std::size_t> modeIndex(std::size_t id) const;  // the position of pair `id`, if present
    Spectrum& setModeIds(std::vector<std::size_t> ids);         // one unique id per pair; empty: 0..n-1
    // Mode `mode` as an Eigenvalue object (continuous rate; multiplier when sampled; conjugate).
    Eigenvalue eigenvalue(std::size_t mode) const;
    Eigenvalue operator[](std::size_t mode) const { return eigenvalue(mode); }
    std::vector<Eigenvalue> eigenvalues() const;

    // Atomic part. Pair `mode`: its eigenvalue, eigenfunction (every representation, each domain bound
    // through `catalog` when it has one), residual and the observables' Koopman modes on it.
    Eigenpair pair(std::size_t mode, const ObservableCatalog* catalog = nullptr) const;
    const std::optional<EigenpairResidual>& residual(std::size_t mode) const;
    bool hasResiduals() const noexcept;
    Spectrum& setResidual(std::size_t mode, std::optional<EigenpairResidual> residual);
    // One residual per mode under one definition and scope (an empty list clears them).
    Spectrum& setResiduals(const std::vector<double>& values, const std::string& definition, const std::string& scope);
    // The eigenfunctions: representations in order (the first is the one a model evaluates), and their
    // values at the training samples (may hold no samples).
    const std::vector<EigenfunctionRepresentation>& representations() const noexcept { return representations_; }
    std::vector<std::string> representationNames() const;
    // By name, else by kind when exactly one representation has it; empty selects the first. An
    // ambiguous kind is refused with the names to choose from.
    const EigenfunctionRepresentation& representation(const std::string& key = {}) const;
    // Adds a representation (a new name, one output per pair), first or last; replaces the one of the
    // same name in place (returning its position); moves one to the front (false when it is there).
    Spectrum& addRepresentation(std::shared_ptr<const FunctionRepresentation> function, bool first = false);
    std::size_t replaceRepresentation(std::shared_ptr<const FunctionRepresentation> function);
    bool moveRepresentationToFront(const std::string& name);
    // The first representation as the same eigenfunctions on an equivalent domain (a bound state map).
    void setFirstRepresentation(std::shared_ptr<const FunctionRepresentation> function);
    const EigenfunctionSamples& samples() const noexcept { return samples_; }
    // Forgets the eigenfunctions (representations and samples), for new ones of the same pairs.
    Spectrum& clearEigenfunctions() { representations_.clear(); samples_ = {}; return *this; }
    Spectrum& setSamples(EigenfunctionSamples samples);
    Eigenfunction eigenfunction(std::size_t mode, const ObservableCatalog* catalog = nullptr) const;

    // Continuous part: bands by name.
    const std::vector<SpectralBand>& bands() const noexcept { return bands_; }
    const SpectralBand* band(const std::string& name) const;
    // Adds a band, or checks that one of the same name is the same set (support, grid, estimate).
    Spectrum& addBand(SpectralBand band);
    Spectrum& setBands(std::vector<SpectralBand> bands);

    // Observables. Koopman modes (outputs x pairs, with offset and fit errors) by observable name.
    const std::map<std::string, ReconstructionResult>& koopmanModes() const noexcept { return koopmanModes_; }
    const ReconstructionResult* koopmanModes(const std::string& observable) const;
    Spectrum& setKoopmanModes(std::string observable, ReconstructionResult modes);
    // Measures, one per (band, observable, output); the band must exist and the density lie
    // on its grid, the packets having one output per grid point.
    const std::vector<SpectralMeasure>& measures() const noexcept { return measures_; }
    const SpectralMeasure* measure(const std::string& band, const std::string& observable, std::size_t output = 0) const;
    Spectrum& addMeasure(SpectralMeasure measure);
    Spectrum& setMeasures(std::vector<SpectralMeasure> measures);
    std::vector<const SpectralMeasure*> bandMeasures(const std::string& band) const;
    // The observable's Koopman modes and measures, as pointers into this spectrum.
    ObservableDecomposition decomposition(const std::string& observable) const;
    // The atoms of the observable's measure, one per pair, by mode id: mass |c_j|^2 ||psi_j||^2 from its
    // Koopman modes, taking ||psi_j|| = 1 (the unit-norm eigenfunctions a spectral fit's normal
    // compression has: eigenvectors of unit length in an orthonormal basis; Colbrook 2023).
    std::vector<SpectralAtom> atoms(const std::string& observable, std::size_t output = 0) const;

    // Evidence and the operator.
    const std::optional<Pseudospectrum>& pseudospectrum() const noexcept { return pseudospectrum_; }
    Spectrum& setPseudospectrum(std::optional<Pseudospectrum> pseudospectrum);
    const KoopmanCompression* compression() const noexcept { return compression_.get(); }
    Spectrum& setCompression(std::shared_ptr<const KoopmanCompression> compression);
    // Functions on spectral sets, by archive key: the pseudospectrum's minimizers
    // (pseudospectrumFunctionKey) and each measure's packets (measureFunctionKey).
    std::vector<std::string> spectralFunctionKeys() const;
    const FunctionExpansion* spectralFunction(const std::string& key) const;
    Spectrum& setSpectralFunction(const std::string& key, std::shared_ptr<const FunctionExpansion> function);
    // Residuals, continuous part, measures and pseudospectrum of `other` (a spectrum of the same size).
    // A measure or pseudospectrum without its functions keeps those this spectrum has for it.
    Spectrum& adoptEvidence(const Spectrum& other);

  private:
    void pairConjugates();
    EigenvalueDomain domain_{EigenvalueDomain::continuousTime};
    double sampleInterval_{};
    std::vector<std::complex<double>> values_, rates_;
    std::vector<bool> rateAvailable_;
    std::vector<int> branches_;
    std::vector<std::optional<std::size_t>> conjugates_;
    std::vector<std::optional<EigenpairResidual>> residuals_;  // empty, or one per mode
    std::vector<std::size_t> modeIds_;                         // empty: id = position
    std::vector<EigenfunctionRepresentation> representations_;
    EigenfunctionSamples samples_;
    std::vector<SpectralBand> bands_;
    std::map<std::string, ReconstructionResult> koopmanModes_;
    std::vector<SpectralMeasure> measures_;
    std::optional<Pseudospectrum> pseudospectrum_;
    std::shared_ptr<const KoopmanCompression> compression_;
    DataOrigin origin_;
};
} // namespace dynsyslab::core
