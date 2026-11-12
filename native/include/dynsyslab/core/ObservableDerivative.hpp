#pragma once
#include <dynsyslab/core/System.hpp>

namespace dynsyslab::core {

// A state map or, when pushForward is true, the final map whose analytic
// Jacobian must be applied to the actual state velocity at each sample.
struct ResolvedObservableDerivative {
    Observable map;
    std::string name, method;
    bool pushForward{};
};

// constantInput is in model units. nullopt denotes a varying input; the
// caller supplies actual velocities when the returned plan requires them.
[[nodiscard]] ResolvedObservableDerivative resolveObservableDerivative(
    const ODESystem& system, const ObservableCatalog& catalog, const std::string& name,
    std::size_t order, std::optional<double> constantInput);
} // namespace dynsyslab::core
