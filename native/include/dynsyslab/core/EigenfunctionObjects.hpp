#pragma once
// Pipeline objects for Koopman spectra and eigenfunctions, independent of the method that found them.
//
//   Spectrum         (core/Spectrum.hpp) the eigenvalues, with the origin that produced them.
//   EigenfunctionSet the eigenfunctions of one spectrum: their values at the training samples, the
//                    origin (method, runtime settings, input data), and an ordered list of
//                    representations. The first representation is the default everywhere,
//                    including error calculations.
//   Representation   a way to evaluate every eigenfunction at a point of its domain. On one
//                    domain it is a FunctionRepresentation (FunctionExpansion or
//                    FunctionInterpolation), which Eigenfunction objects point at.
//                    SharedRepresentation wraps any one of them; interpolations are built by the
//                    same interpolateSamples every Observable uses.
#include <dynsyslab/core/Spectrum.hpp>
#include <dynsyslab/core/Eigenfunction.hpp>
#include <dynsyslab/core/Eigenvalue.hpp>
#include <dynsyslab/core/FunctionInterpolation.hpp>
#include <dynsyslab/core/FunctionRepresentation.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/TimeSeriesRecord.hpp>

#include <complex>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace dynsyslab::core {

// Values of every eigenfunction at the training samples the method produced them on.
struct EigenfunctionSamples {
    std::vector<double> time;               // label time of each sample
    std::vector<std::size_t> sampleIndices;  // row of the input record each sample is labelled by
    std::vector<std::complex<double>> values;  // samples x modes, row-major
    std::size_t modes{};
    [[nodiscard]] std::size_t samples() const noexcept { return time.size(); }
};

class EigenfunctionRepresentation {
public:
    virtual ~EigenfunctionRepresentation() = default;
    [[nodiscard]] virtual std::string kind() const = 0;
    [[nodiscard]] virtual std::vector<std::string> domains() const = 0;  // first is the default
    [[nodiscard]] virtual std::size_t domainDimension(const std::string& domain) const = 0;
    // Identity of the domain's coordinates, matched against TimeSeriesRecord channel observables.
    [[nodiscard]] virtual std::vector<std::string> domainObservables(const std::string& domain) const = 0;
    [[nodiscard]] virtual std::vector<std::complex<double>> evaluate(const double* point,
                                                                     const std::string& domain) const = 0;
    // Representation values at every training sample (samples x modes). `input` is the record the
    // eigenfunctions were fitted on, for representations that rebuild their argument from it.
    [[nodiscard]] virtual std::vector<std::complex<double>> evaluateTraining(
        const EigenfunctionSamples& samples, const std::string& domain, const TimeSeriesRecord* input) const = 0;
    [[nodiscard]] virtual std::shared_ptr<const EigenfunctionRepresentation> extract(
        const std::vector<std::size_t>& modes) const = 0;
    [[nodiscard]] virtual json::Value describe() const = 0;
    // The representation on one domain, for every mode (no eigenvalues).
    [[nodiscard]] virtual std::shared_ptr<const FunctionRepresentation> function(const std::string& domain) const = 0;
};

// Every eigenfunction on one domain through one FunctionRepresentation: the general representation,
// independent of the method that produced the eigenfunctions. An interpolation built by
// interpolateSamples (FunctionInterpolation.hpp) carries its domain and base points; any method
// (H-DMD, EDMD, SINDy, ...) uses it the same way an Observable does.
class SharedRepresentation final : public EigenfunctionRepresentation {
public:
    explicit SharedRepresentation(std::shared_ptr<const FunctionRepresentation> function);
    [[nodiscard]] std::string kind() const override { return function_->kind(); }
    // One domain, named "state" for the full state and after its coordinates otherwise.
    [[nodiscard]] std::vector<std::string> domains() const override { return {domain_}; }
    [[nodiscard]] std::size_t domainDimension(const std::string& domain) const override;
    [[nodiscard]] std::vector<std::string> domainObservables(const std::string& domain) const override;
    [[nodiscard]] std::vector<std::complex<double>> evaluate(const double* point, const std::string& domain) const override;
    // An interpolation evaluates at its stored base points; other functions cannot rebuild them.
    [[nodiscard]] std::vector<std::complex<double>> evaluateTraining(
        const EigenfunctionSamples& samples, const std::string& domain, const TimeSeriesRecord* input) const override;
    [[nodiscard]] std::shared_ptr<const EigenfunctionRepresentation> extract(const std::vector<std::size_t>& modes) const override;
    [[nodiscard]] json::Value describe() const override;
    [[nodiscard]] std::shared_ptr<const FunctionRepresentation> function(const std::string& domain) const override;
    [[nodiscard]] const std::shared_ptr<const FunctionRepresentation>& function() const noexcept { return function_; }

private:
    void requireDomain(const std::string& domain) const;
    std::shared_ptr<const FunctionRepresentation> function_;
    std::string domain_;
};

// The name SharedRepresentation gives a domain: "state" for the full state, else its coordinate
// names joined by commas, else "coordinates".
[[nodiscard]] std::string domainName(const FunctionDomain& domain);

struct EigenfunctionErrors {
    std::string representation, domain;
    std::size_t samples{};
    std::vector<double> valueRelative;     // ||R(p_k) - psi_k|| / ||psi_k|| over the samples, per mode
    std::vector<double> koopmanResidual;   // ||R(p_k+1) - mu^dt R(p_k)|| / ||R(p_k)||, per mode
};

class EigenfunctionSet {
public:
    EigenfunctionSet(Spectrum spectrum, EigenfunctionSamples samples, DataOrigin origin,
                     std::vector<std::shared_ptr<const EigenfunctionRepresentation>> representations);
    [[nodiscard]] const Spectrum& spectrum() const noexcept { return spectrum_; }
    [[nodiscard]] const EigenfunctionSamples& samples() const noexcept { return samples_; }
    [[nodiscard]] const DataOrigin& origin() const noexcept { return origin_; }
    [[nodiscard]] std::size_t size() const noexcept { return spectrum_.size(); }
    [[nodiscard]] const std::vector<std::shared_ptr<const EigenfunctionRepresentation>>& representations() const noexcept {
        return representations_;
    }
    // By kind; empty selects the first (default) representation.
    [[nodiscard]] const EigenfunctionRepresentation& representation(const std::string& kind = {}) const;
    // The domain of `representation` whose observables match `observables`, else its default domain
    // when `observables` is empty. Throws when nothing matches.
    [[nodiscard]] std::string recognizeDomain(const std::vector<std::string>& observables,
                                              const std::string& representation = {}) const;
    [[nodiscard]] std::vector<std::complex<double>> evaluate(const double* point, const std::string& domain = {},
                                                             const std::string& representation = {}) const;
    [[nodiscard]] EigenfunctionErrors errors(const std::string& representation = {}, const std::string& domain = {},
                                             const TimeSeriesRecord* input = nullptr) const;
    [[nodiscard]] EigenfunctionSet extract(const std::vector<std::size_t>& modes) const;
    // Mode `mode` as an Eigenfunction object: the chosen representation on the chosen domain
    // (defaults: the first ones) and the spectrum's eigenvalue, through the domain's state map
    // (`catalog` binds a catalog map); a function of q when the domain has none.
    [[nodiscard]] Eigenfunction eigenfunction(std::size_t mode, const std::string& domain = {},
                                              const std::string& representation = {},
                                              const ObservableCatalog* catalog = nullptr) const;
    // Replaces the representation of the same kind, keeping its position in the order.
    void replaceRepresentation(std::shared_ptr<const EigenfunctionRepresentation> representation);

private:
    Spectrum spectrum_;
    EigenfunctionSamples samples_;
    DataOrigin origin_;
    std::vector<std::shared_ptr<const EigenfunctionRepresentation>> representations_;
};

// The set's interpolation, and switching its interpolant in place (for sweeps: set the method,
// evaluate errors, set another). The domain and base points stay.
[[nodiscard]] const FunctionInterpolation& interpolationOf(const EigenfunctionSet& set);
void setInterpolationMethod(EigenfunctionSet& set, InterpolationMethod method, std::size_t neighbors);

}  // namespace dynsyslab::core
