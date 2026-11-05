#pragma once

#include <dynsyslab/core/AnalyticExpression.hpp>
#include <complex>
#include <cstddef>
#include <functional>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace dynsyslab::core {

enum class ObservableJacobianBackend { automatic, analytic, finiteDifference };

struct ObservableDescriptor {
    std::string name;
    std::string label;
    std::string units;
    std::string stateCoordinateSystem;
    std::vector<std::string> outputNames;
    bool complexValued{};
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
    [[nodiscard]] std::size_t stateDimension() const noexcept;
    [[nodiscard]] std::size_t outputDimension() const noexcept;
    // Exact y[output] = state[index] provenance, assigned only by the defining factories.
    // Custom evaluators remain unknown, regardless of names or sampled values.
    [[nodiscard]] std::optional<std::size_t> stateComponentIndex(std::size_t output) const;
    [[nodiscard]] bool hasAnalyticJacobian() const noexcept;
    [[nodiscard]] const std::vector<AnalyticExpression>& expressions() const noexcept { return expressions_; }
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
    std::vector<AnalyticExpression> expressions_;
    Evaluator evaluator_;
    Jacobian analyticJacobian_;
    std::vector<std::optional<std::size_t>> stateComponents_;
};

class ObservableCatalog {
public:
    ObservableCatalog(std::vector<Observable> observables,
        std::map<std::string, std::string> aliases = {},
        std::map<std::string, std::vector<std::string>> coordinatePresets = {});

    [[nodiscard]] bool contains(const std::string& name) const;
    [[nodiscard]] const Observable& at(const std::string& name) const;
    [[nodiscard]] Observable coordinates(const std::vector<std::string>& names,
                                          std::string name = {}) const;
    [[nodiscard]] Observable coordinatePreset(const std::string& name) const;
    [[nodiscard]] std::vector<std::string> names() const;
    [[nodiscard]] const std::map<std::string, std::string>& aliases() const noexcept;
    [[nodiscard]] const std::map<std::string, std::vector<std::string>>&
        coordinatePresets() const noexcept;

private:
    std::map<std::string, Observable> observables_;
    std::map<std::string, std::string> aliases_;
    std::map<std::string, std::vector<std::string>> coordinatePresets_;
};

}  // namespace dynsyslab::core
