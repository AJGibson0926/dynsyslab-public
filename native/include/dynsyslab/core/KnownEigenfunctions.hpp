#pragma once
#include <dynsyslab/core/PolynomialEigenfunctions.hpp>
#include <dynsyslab/core/System.hpp>
#include <optional>

namespace dynsyslab::core {
// Exact, system-specific mathematics; never a fitted bank or sampled certificate.
// Owned by the eigenfunction module, so systems do not depend on Koopman runtime.
struct KnownEigenfunctionResult {
    std::string identifier, status{"unavailable"}, reason, conditions;
    std::vector<std::string> names, references;
    std::optional<PolynomialEigenfunctionSet> eigenfunctions;
    // coordinate x output mode x predictor mode; partial_j z = D_j z.
    std::vector<double> derivativeOperators;
    // physical state x mode, if a linear state readout exists.
    std::vector<double> stateReadout;
};
KnownEigenfunctionResult knownPolynomialEigenfunctions(const ODESystem& system);
json::Value knownEigenfunctionsJson(const KnownEigenfunctionResult& result);
} // namespace dynsyslab::core
