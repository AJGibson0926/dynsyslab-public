#pragma once
// The spectrum of a set of Koopman eigenfunctions, with the origin (method, settings, inputs) that
// produced it. It has two parts and optional evidence:
//   point part       the eigenvalues, one per eigenpair, each with the pair's residual when the method
//                    computed one (Eigenpair.hpp). The eigenfunctions themselves are held by the
//                    KoopmanModel that shares their representations; KoopmanModel::eigenpair(j) joins
//                    them, `pair(j)` here gives the eigenvalue and residual.
//   continuous part  components with a spectral measure or density and no eigenfunctions
//                    (SpectralMeasure.hpp); empty unless the method computed one (mpEDMD, Rigged DMD,
//                    ResDMD measures for measure-preserving systems).
//   pseudospectrum   optional computed evidence: the minimum residual on a grid of points z
//                    (Pseudospectrum.hpp; Trefethen & Embree 2005; Colbrook & Townsend 2024).
// Selecting modes keeps their residuals; the continuous part and pseudospectrum describe the
// operator, not the selection, and are kept as well.
//
// Each mode stores its value in its own time domain (a continuous rate lambda, or a sampled
// multiplier mu = exp(lambda dt)) plus the selected logarithm branch, so conversion is explicit and
// lossless. A zero sampled multiplier is valid but has no finite continuous rate. `eigenvalue(j)`
// and `eigenvalues()` hand out Eigenvalue objects, with conjugate partners paired.
#include <dynsyslab/core/Eigenpair.hpp>
#include <dynsyslab/core/Pseudospectrum.hpp>
#include <dynsyslab/core/SpectralMeasure.hpp>
#include <dynsyslab/core/TimeSeriesRecord.hpp>

#include <complex>
#include <cstddef>
#include <optional>
#include <vector>

namespace dynsyslab::core {
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
    Spectrum extract(const std::vector<std::size_t>& modes) const;
    // Mode `mode` as an Eigenvalue object (continuous rate; multiplier when sampled; conjugate).
    Eigenvalue eigenvalue(std::size_t mode) const;
    Eigenvalue operator[](std::size_t mode) const { return eigenvalue(mode); }
    std::vector<Eigenvalue> eigenvalues() const;

    // Point part: mode `mode`'s eigenvalue and residual (no eigenfunction; see KoopmanModel::eigenpair).
    Eigenpair pair(std::size_t mode) const;
    const std::optional<EigenpairResidual>& residual(std::size_t mode) const;
    bool hasResiduals() const noexcept;
    Spectrum& setResidual(std::size_t mode, std::optional<EigenpairResidual> residual);
    // One residual per mode under one definition and scope (an empty list clears them).
    Spectrum& setResiduals(const std::vector<double>& values, const std::string& definition, const std::string& scope);
    // Continuous part and pseudospectrum.
    const std::vector<SpectralMeasure>& continuous() const noexcept { return continuous_; }
    Spectrum& setContinuous(std::vector<SpectralMeasure> components) { continuous_ = std::move(components); return *this; }
    const std::optional<Pseudospectrum>& pseudospectrum() const noexcept { return pseudospectrum_; }
    Spectrum& setPseudospectrum(std::optional<Pseudospectrum> pseudospectrum);
    // Residuals, continuous part and pseudospectrum of `other` (a spectrum of the same size).
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
    std::vector<SpectralMeasure> continuous_;
    std::optional<Pseudospectrum> pseudospectrum_;
    DataOrigin origin_;
};
} // namespace dynsyslab::core
