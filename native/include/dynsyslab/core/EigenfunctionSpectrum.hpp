#pragma once
#include <complex>
#include <cstddef>
#include <optional>
#include <vector>

namespace dynsyslab::core {
enum class EigenvalueDomain { continuousTime, discreteTime };

// Stores the asserted spectrum in its own time domain. Conversion is explicit;
// zero sampled multipliers are valid but have no finite continuous generator.
class EigenfunctionSpectrum {
  public:
    static EigenfunctionSpectrum fromContinuousRates(std::vector<std::complex<double>> rates,
                                                     double sampleInterval = 0);
    static EigenfunctionSpectrum
    fromDiscreteMultipliers(std::vector<std::complex<double>> multipliers, double sampleInterval,
                            std::vector<int> logarithmBranches = {});
    // Compatibility import: preserve archived rates exactly, including their
    // selected branches, rather than round-tripping through exp/log.
    static EigenfunctionSpectrum fromStoredRates(std::vector<std::complex<double>> rates,
                                                 double sampleInterval, EigenvalueDomain origin,
                                                 std::vector<int> logarithmBranches = {});
    // Lossless archive round-trip: validate and retain both recorded values.
    static EigenfunctionSpectrum fromStoredDiscrete(std::vector<std::complex<double>> multipliers,
                                                    double sampleInterval,
                                                    std::vector<int> logarithmBranches,
                                                    const std::vector<std::complex<double>> &rates);
    EigenvalueDomain domain() const noexcept { return domain_; }
    double sampleInterval() const noexcept { return sampleInterval_; }
    std::size_t size() const noexcept { return values_.size(); }
    const std::vector<std::complex<double>> &values() const noexcept { return values_; }
    const std::vector<int> &logarithmBranches() const noexcept { return branches_; }
    bool hasContinuousRates() const noexcept;
    std::optional<std::complex<double>> continuousRate(std::size_t mode) const;
    const std::vector<std::complex<double>> &continuousRates() const;
    std::vector<std::complex<double>> discreteMultipliers(double sampleInterval = 0) const;
    // Nonzero sampled modes use the recorded logarithm branch between samples.
    // Zero multipliers only admit nonnegative integer sampling intervals.
    std::complex<double> transition(std::size_t mode, double elapsed) const;
    EigenfunctionSpectrum extract(const std::vector<std::size_t> &modes) const;

  private:
    EigenfunctionSpectrum() = default;
    EigenvalueDomain domain_{EigenvalueDomain::continuousTime};
    double sampleInterval_{};
    std::vector<std::complex<double>> values_, rates_;
    std::vector<bool> rateAvailable_;
    std::vector<int> branches_;
};
} // namespace dynsyslab::core
