#pragma once

#include <dynsyslab/core/TimeSeries.hpp>
#include <dynsyslab/core/AnalyticExpression.hpp>
#include <complex>
#include <cstddef>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace dynsyslab::core {

class FunctionRepresentation;
class ObservableCatalog;

enum class ObservableJacobianBackend { automatic, analytic, finiteDifference };

// Passive provenance only. The Koopman/workflow layers resolve and validate a reference.
// An artifact locator is not an identity: the frozen model SHA-256 is required.
struct ObservableDecompositionReference {
    std::string artifact, modelIdentity, component;
    std::string kind{"reconstruction"}; // reconstruction, eigenfunction, spectral-measure, wave-packet
    std::size_t index{};
};
struct ObservableDescriptor {
    std::string name;
    std::string label;
    std::string units;
    std::string stateCoordinateSystem;
    std::vector<std::string> outputNames;
    bool complexValued{};
    std::vector<ObservableDecompositionReference> decompositions;
};

class Observable {
public:
    using Values = std::vector<std::complex<double>>;
    using Evaluator = std::function<Values(const double*)>;
    using Jacobian = std::function<Values(const double*)>;

    Observable(ObservableDescriptor descriptor, std::size_t stateDimension,
               std::size_t outputDimension, Evaluator evaluator,
               Jacobian analyticJacobian = {});

    [[nodiscard]] const ObservableDescriptor& descriptor() const noexcept;
    [[nodiscard]] Observable withDecomposition(ObservableDecompositionReference reference) const;
    [[nodiscard]] std::size_t stateDimension() const noexcept;
    [[nodiscard]] std::size_t outputDimension() const noexcept;
    // Exact y[output] = state[index] provenance, assigned only by the defining factories.
    // Custom evaluators remain unknown, regardless of names or sampled values.
    [[nodiscard]] std::optional<std::size_t> stateComponentIndex(std::size_t output) const;
    [[nodiscard]] bool hasAnalyticJacobian() const noexcept;
    // The closed form, one real expression per output: from the expansion representation (on the
    // physical state, every term with a closed form) or from state components; empty otherwise.
    [[nodiscard]] std::vector<AnalyticExpression> expressions() const;
    [[nodiscard]] Observable partialDerivative(std::size_t axis) const;
    [[nodiscard]] Observable flowDerivative(const std::vector<AnalyticExpression>& field,
                                           std::size_t order = 1) const;
    // Total derivatives with supplied affine input; rates are needed at order two.
    // Higher orders require held-constant inputs (zero inputRates).
    [[nodiscard]] Observable inputFlowDerivative(
        const std::vector<AnalyticExpression>& drift,
        const std::vector<AnalyticExpression>& inputFields,
        const std::vector<double>& inputs, const std::vector<double>& inputRates = {},
        std::size_t order = 1) const;
    [[nodiscard]] static Observable fromExpressions(ObservableDescriptor descriptor,
        std::vector<AnalyticExpression> expressions);
    [[nodiscard]] json::Value expressionDefinition() const;
    [[nodiscard]] static Observable fromExpressionDefinition(const json::Value& value);
    [[nodiscard]] Values evaluate(const double* state) const;
    // The observable along a time series, as another time series: at each sample's state, or, when
    // its domain is a series window (FunctionDomain::window), at the window cut at each sample;
    // samples without a complete window are skipped. A complex-valued observable gives a real and an
    // imaginary column per output. Eigenfunctions inherit this.
    [[nodiscard]] TimeSeries evaluate(const TimeSeries& series) const;
    [[nodiscard]] Values jacobian(const double* state,
        ObservableJacobianBackend backend = ObservableJacobianBackend::automatic,
        double relativeStep = 1.0e-6) const;
    [[nodiscard]] std::vector<double> evaluateReal(const double* state) const;
    [[nodiscard]] std::vector<double> evaluateRealBatch(
        const std::vector<double>& states, std::size_t samples) const;

    // Push supplied state velocities through the observable Jacobian, sample by sample.
    [[nodiscard]] std::vector<double> pushForwardRealBatch(
        const std::vector<double>& states, const std::vector<double>& velocities,
        std::size_t samples,
        ObservableJacobianBackend backend = ObservableJacobianBackend::automatic) const;

    // g(q(x)) from a representation (expansion or interpolation) on its domain q, through the
    // domain's state map q(x); a function of q itself when the domain has none. `catalog` binds a
    // catalog map that is known by name only. The Jacobian is the exact chain rule when the
    // representation and the map are exact.
    [[nodiscard]] static Observable fromRepresentation(ObservableDescriptor descriptor,
        std::shared_ptr<const FunctionRepresentation> representation,
        const ObservableCatalog* catalog = nullptr);
    // The representation property: an expansion (formulas included: fromExpressions writes each
    // formula as a dictionary term) or an interpolation; null for evaluator observables.
    [[nodiscard]] const std::shared_ptr<const FunctionRepresentation>& representation() const noexcept { return representation_; }
    // "expansion", "interpolation" or "evaluator".
    [[nodiscard]] std::string representationKind() const;

    [[nodiscard]] static Observable stateComponent(
        std::size_t stateDimension, std::size_t index,
        std::string name, std::string label,
        std::string coordinates = "native");
    [[nodiscard]] static Observable identity(
        std::size_t stateDimension, std::vector<std::string> names,
        std::vector<std::string> labels, std::string coordinates = "native");
    [[nodiscard]] static Observable stack(
        std::string name, std::string label,
        const std::vector<Observable>& observables);

private:
    ObservableDescriptor descriptor_;
    std::size_t stateDimension_{};
    std::size_t outputDimension_{};
    Evaluator evaluator_;
    Jacobian analyticJacobian_;
    std::vector<std::optional<std::size_t>> stateComponents_;
    std::shared_ptr<const FunctionRepresentation> representation_;
};

// Named first derivative on the system's native clock. Some formulas are
// Lie derivatives along the zero-input drift only.
struct ObservableTimeDerivative {
    std::string target;
    bool zeroInputOnly{};
};
using ObservableDerivativeLinks = std::map<std::string, ObservableTimeDerivative>;

class ObservableCatalog {
public:
    ObservableCatalog(std::vector<Observable> observables,
        std::map<std::string, std::string> aliases = {},
        std::map<std::string, std::vector<std::string>> coordinatePresets = {},
        ObservableDerivativeLinks derivatives = {});

    [[nodiscard]] bool contains(const std::string& name) const;
    [[nodiscard]] const Observable& at(const std::string& name) const;
    [[nodiscard]] const ObservableTimeDerivative* timeDerivative(const std::string& name) const;
    [[nodiscard]] const ObservableDerivativeLinks& timeDerivatives() const noexcept { return derivatives_; }
    [[nodiscard]] Observable coordinates(const std::vector<std::string>& names,
                                          std::string name = {}) const;
    [[nodiscard]] Observable coordinatePreset(const std::string& name) const;
    [[nodiscard]] std::vector<std::string> names() const;
    [[nodiscard]] const std::map<std::string, std::string>& aliases() const noexcept;
    [[nodiscard]] const std::map<std::string, std::vector<std::string>>&
        coordinatePresets() const noexcept;

private:
    ObservableDerivativeLinks derivatives_;
    std::map<std::string, Observable> observables_;
    std::map<std::string, std::string> aliases_;
    std::map<std::string, std::vector<std::string>> coordinatePresets_;
};

}  // namespace dynsyslab::core
