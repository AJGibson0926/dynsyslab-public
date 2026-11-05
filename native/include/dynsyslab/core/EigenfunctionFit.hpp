#pragma once
#include <dynsyslab/core/PolynomialEigenfunctions.hpp>
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
struct EigenfunctionFitResult {
    PolynomialEigenfunctionSet eigenfunctions;
    EigenfunctionFitDiagnostics diagnostics;
};
// Explicit compatibility boundary for the historical algorithm-named DTO.
PolynomialEigenfunctionSet polynomialEigenfunctionsFromLegacy(const KoopmanSINDyResult &);
KoopmanSINDyResult legacyPolynomialEigenfunctions(const PolynomialEigenfunctionSet &,
                                                  const EigenfunctionFitDiagnostics & = {});
EigenfunctionFitDiagnostics fitDiagnosticsFromLegacy(const KoopmanSINDyResult &);

// Same generator fit as fitKoopmanEigenfunctions; coordinates are explicit and
// diagnostics are separate. Supplied derivatives may be measured or modeled.
EigenfunctionFitResult fitSindyEigenfunctions(
    const std::vector<double> &values, std::size_t samples, const std::vector<double> &time,
    const std::vector<std::complex<double>> &continuousRates, EigenfunctionCoordinates coordinates,
    const std::vector<double> &derivatives = {}, const KoopmanSINDyOptions &options = {});
} // namespace dynsyslab::core
