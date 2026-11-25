#pragma once
#include <dynsyslab/core/Spectrum.hpp>
#include <dynsyslab/core/FunctionExpansion.hpp>
#include <dynsyslab/core/SINDy.hpp>

namespace dynsyslab::core {
// Evidence from a fitting procedure, NOT a property of the polynomial itself.
// Empty algorithm/metric means unspecified, not a claim of zero error.
struct EigenfunctionFitDiagnostics {
    std::string algorithm, residualDefinition;
    std::size_t sampleCount{}, directFitCount{};
    std::vector<double> residualRms, residualMaximum;
    EigenfunctionFitDiagnostics extract(const std::vector<std::size_t> &modes) const;
};
// What any producer (a fit, a decomposition, a library) hands back: the eigenfunctions' shared
// expansion and, separately, their eigenvalues. KoopmanModel and Eigenfunction objects are built
// from the pair; the expansion itself carries no eigenvalues.
struct EigenfunctionFitResult {
    FunctionExpansion functions;
    Spectrum spectrum;
    EigenfunctionFitDiagnostics diagnostics;
};
// Explicit compatibility boundary for the historical algorithm-named DTO.
EigenfunctionFitResult polynomialEigenfunctionsFromLegacy(const KoopmanSINDyResult &);
KoopmanSINDyResult legacyPolynomialEigenfunctions(const FunctionExpansion &, const Spectrum &,
                                                  const EigenfunctionFitDiagnostics & = {});
EigenfunctionFitDiagnostics fitDiagnosticsFromLegacy(const KoopmanSINDyResult &);

// Same generator fit as fitKoopmanEigenfunctions; coordinates are explicit and
// diagnostics are separate. Supplied derivatives may be measured or modeled.
EigenfunctionFitResult fitSindyEigenfunctions(
    const std::vector<double> &values, std::size_t samples, const std::vector<double> &time,
    const std::vector<std::complex<double>> &continuousRates, FunctionDomain coordinates,
    const std::vector<double> &derivatives = {}, const KoopmanSINDyOptions &options = {});
} // namespace dynsyslab::core
