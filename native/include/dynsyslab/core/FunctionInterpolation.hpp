#pragma once
// The interpolation representation: outputs interpolated from sampled values at scattered points of
// a domain, with a selectable interpolant (InterpolationMethod). Like FunctionExpansion it is a
// FunctionRepresentation and holds no eigenvalues; Observables and Eigenfunctions point at it.
#include <dynsyslab/core/FunctionRepresentation.hpp>

#include <complex>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

#include <dynsyslab/core/Json.hpp>

namespace dynsyslab::core {

// An interpolant method for scattered samples (see ScatteredInterpolants.cpp). Every method works on
// the k nearest samples; `neighbors` lives with the data.
struct InterpolationMethod {
    std::string name{"inverse_distance"};
    double power{2};   // inverse_distance: weights 1/d^power
    double width{0};   // gaussian_kernel / local_linear bandwidth, RBF shape length; 0 = mean neighbour distance
};
[[nodiscard]] std::vector<std::string> interpolationMethods();
[[nodiscard]] json::Value interpolationMethodJson(const InterpolationMethod& method, std::size_t neighbors);
// Reads a method object; writes its "neighbors" (when given) to *neighbors.
[[nodiscard]] InterpolationMethod interpolationMethodFromJson(const json::Value& value, std::size_t* neighbors = nullptr);

struct InterpolationData {
    std::vector<double> trainingStates;                 // sample x domain coordinate
    // Empty: interpolate on the physical state. Otherwise the catalog observable
    // whose values are the interpolation coordinates, e.g. psi(echo).
    std::string coordinateObservable;
    std::vector<std::complex<double>> trainingValues;  // sample x output
    std::size_t samples{};
    std::size_t stateDimension{};
    std::size_t neighbors{8};
    InterpolationMethod method;
    double finiteDifferenceRelativeStep{1.0e-6};
    std::string provenance;
    [[nodiscard]] std::size_t outputs() const noexcept { return samples ? trainingValues.size() / samples : 0; }
};

// Immutable multi-output interpolant. One spatial search is shared by every output, which keeps
// large H-DMD K-ROM evaluation much cheaper than one interpolator/tree per eigenfunction.
class FunctionInterpolation final : public FunctionRepresentation {
public:
    explicit FunctionInterpolation(InterpolationData data);
    // On an explicit domain (its dimension must match the samples').
    FunctionInterpolation(InterpolationData data, FunctionDomain coordinates);

    [[nodiscard]] std::string kind() const override { return "interpolation"; }
    [[nodiscard]] const FunctionDomain& coordinates() const noexcept override { return coordinates_; }
    [[nodiscard]] std::size_t outputCount() const noexcept override { return data_.outputs(); }
    [[nodiscard]] const InterpolationData& data() const noexcept;
    [[nodiscard]] std::size_t modes() const noexcept { return outputCount(); }
    [[nodiscard]] std::vector<std::complex<double>> evaluate(const double* point) const override;
    // Row-major output x coordinate gradient (central differences).
    [[nodiscard]] std::vector<std::complex<double>> jacobian(const double* point) const override;
    [[nodiscard]] std::vector<std::complex<double>> gradient(const double* point) const { return jacobian(point); }
    [[nodiscard]] bool exactJacobian() const noexcept override { return false; }
    [[nodiscard]] std::shared_ptr<const FunctionInterpolation> select(const std::vector<std::size_t>& outputs) const;
    [[nodiscard]] std::shared_ptr<const FunctionRepresentation> selectOutputs(
        const std::vector<std::size_t>& outputs) const override { return select(outputs); }
    [[nodiscard]] std::shared_ptr<const FunctionRepresentation> onDomain(FunctionDomain coordinates) const override;
    // The same samples and base points with another interpolant.
    [[nodiscard]] std::shared_ptr<const FunctionInterpolation> withMethod(InterpolationMethod method,
                                                                         std::size_t neighbors) const;
    // The base points q_k (samples x dimension), stored in the interpolation.
    [[nodiscard]] const std::vector<double>& points() const noexcept { return data_.trainingStates; }
    [[nodiscard]] std::shared_ptr<const FunctionInterpolation> extract(const std::vector<std::size_t>& outputs) const {
        return select(outputs);
    }

private:
    struct Index;
    InterpolationData data_;
    FunctionDomain coordinates_;
    std::shared_ptr<const Index> index_;
};

// Sampled values interpolated on a domain q = q(x): the one way any Observable, and so any
// Eigenfunction, is built from samples. `states` are the full states x_k at the samples
// (samples x physical dimension), which the caller finds through the sampled series' origin; the
// domain's state map turns them into the base points q_k = q(x_k), which the interpolation stores.
// The default domain is the full state (identity map). A domain without a state map takes `states`
// as the base points themselves (measured or history coordinates). A catalog map must be bound
// (FunctionDomain::bound) first. `values` are samples x outputs. Each observable chooses its own
// interpolant: `method` and `neighbors` belong to the result, not to the framework.
[[nodiscard]] std::shared_ptr<const FunctionInterpolation> interpolateSamples(
    const std::vector<double>& states, std::size_t samples, const FunctionDomain& domain,
    std::vector<std::complex<double>> values, InterpolationMethod method = {}, std::size_t neighbors = 8,
    std::string provenance = {});

// The configured interpolant at `point` from its nearest samples (squared distance, sample index).
[[nodiscard]] std::vector<std::complex<double>> interpolateScatteredNeighbours(const InterpolationData& data,
    const std::vector<std::pair<double, std::size_t>>& neighbours, const double* point);

}  // namespace dynsyslab::core
