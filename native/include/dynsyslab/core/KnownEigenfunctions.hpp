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
//   libraryEigenfunctions      one Eigenfunction object per entry with a constant eigenvalue, on the
//                              physical state
//   libraryEigenfunctionBank   the same as one shared FunctionExpansion (outputs named by entry id,
//                              terms from FunctionExpansion::fromExpressions) plus their spectrum,
//                              whose bands hold the shell entries requested (all, by default)
//   libraryBands               the "shell" entries as spectral bands
// A "shell" entry e^{i n phi} has eigenvalue i n Omega(H), constant on each level set of an invariant H
// but not globally: where the frequency varies (twist) it is not an eigenfunction on the phase space
// but the profile of generalized eigenfunctions delta(n Omega(H(x)) - p) e^{i n phi(x)} on a band of the
// imaginary axis (Mezic 2020, section 11), parametrized by the frequency p, with level map n Omega(H).
// `grid` samples the band's frequencies when given.
[[nodiscard]] std::vector<Eigenfunction> libraryEigenfunctions(const ODESystem& system,
                                                               const std::vector<std::string>& ids = {});
[[nodiscard]] EigenfunctionFitResult libraryEigenfunctionBank(const ODESystem& system,
                                                              const std::vector<std::string>& ids = {});
[[nodiscard]] std::vector<SpectralBand> libraryBands(const ODESystem& system, const std::vector<std::string>& ids = {},
                                                     const std::vector<double>& grid = {});
// Per eigenfunction, max |grad psi . f - lambda psi| / max |psi| over `states` (rows of the physical
// state), f the unforced drift: zero for an exact eigenfunction, small near x* for a linearized
// principal one. Works for any Eigenfunction (expansion or interpolation); checks library entries
// against the configured equations.
[[nodiscard]] std::vector<double> eigenfunctionGeneratorResidual(const ODESystem& system,
                                                                 const std::vector<Eigenfunction>& eigenfunctions,
                                                                 const std::vector<double>& states);
json::Value knownEigenfunctionsJson(const KnownEigenfunctionResult& result);
} // namespace dynsyslab::core
