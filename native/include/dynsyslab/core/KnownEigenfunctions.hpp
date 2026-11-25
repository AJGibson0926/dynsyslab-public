#pragma once
#include <dynsyslab/core/Eigenfunction.hpp>
#include <dynsyslab/core/EigenfunctionFit.hpp>
#include <dynsyslab/core/System.hpp>
#include <optional>

namespace dynsyslab::core {
// Exact, system-specific mathematics; never a fitted bank or sampled certificate.
// Owned by the eigenfunction module, so systems do not depend on Koopman runtime.
struct KnownEigenfunctionResult {
    std::string identifier, status{"unavailable"}, reason, conditions;
    std::vector<std::string> names, references;
    std::optional<EigenfunctionFitResult> eigenfunctions;  // shared expansion + spectrum
    // coordinate x output mode x predictor mode; partial_j z = D_j z.
    std::vector<double> derivativeOperators;
    // physical state x mode, if a linear state readout exists.
    std::vector<double> stateReadout;
};
KnownEigenfunctionResult knownPolynomialEigenfunctions(const ODESystem& system);

// The system's own eigenfunction library (DynamicalKnowledge::eigenfunctions, hand-derived exact
// and principal eigenfunctions) with the configured parameters substituted. `ids` selects entries
// in that order; by default every applicable entry. Throws if a requested entry is missing or not
// applicable.
//   libraryEigenfunctions      one Eigenfunction object per entry, on the physical state
//   libraryEigenfunctionBank   the same as one shared FunctionExpansion (outputs named by entry id,
//                              terms from FunctionExpansion::fromExpressions) plus their spectrum,
//                              e.g. to build a KoopmanModel or project observables
// A "shell" entry has an eigenvalue that is constant on each level set of an invariant (e.g.
// i n Omega(H) for a phase e^{i n phi}); it is an eigenfunction along each orbit, and the eigenvalue
// is taken on the orbit through `reference` (a physical state), which shell entries require.
[[nodiscard]] std::vector<Eigenfunction> libraryEigenfunctions(const ODESystem& system,
                                                               const std::vector<std::string>& ids = {},
                                                               const std::vector<double>& reference = {});
[[nodiscard]] EigenfunctionFitResult libraryEigenfunctionBank(const ODESystem& system,
                                                              const std::vector<std::string>& ids = {},
                                                              const std::vector<double>& reference = {});
// Per eigenfunction, max |grad psi . f - lambda psi| / max |psi| over `states` (rows of the physical
// state), f the unforced drift: zero for an exact eigenfunction, small near x* for a linearized
// principal one. Works for any Eigenfunction (expansion or interpolation); checks library entries
// against the configured equations.
[[nodiscard]] std::vector<double> eigenfunctionGeneratorResidual(const ODESystem& system,
                                                                 const std::vector<Eigenfunction>& eigenfunctions,
                                                                 const std::vector<double>& states);
json::Value knownEigenfunctionsJson(const KnownEigenfunctionResult& result);
} // namespace dynsyslab::core
