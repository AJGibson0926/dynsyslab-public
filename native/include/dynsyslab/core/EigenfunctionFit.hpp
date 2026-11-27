#pragma once
#include <dynsyslab/core/Spectrum.hpp>
#include <dynsyslab/core/FunctionExpansion.hpp>

#include <optional>

namespace dynsyslab::core {
// Evidence from a fitting procedure, NOT a property of the polynomial itself.
// Empty algorithm/metric means unspecified, not a claim of zero error.
// The pairs' residuals themselves are on the result's Spectrum (EigenpairResidual, Eigenpair.hpp); the
// fit keeps only the largest pointwise residual under the same definition.
struct EigenfunctionFitDiagnostics {
    std::string algorithm;
    std::size_t sampleCount{}, directFitCount{};
    std::vector<double> residualMaximum;
    EigenfunctionFitDiagnostics extract(const std::vector<std::size_t> &modes) const;
};
// What any producer (a fit, a decomposition, a library) hands back: the eigenfunctions' shared
// expansion and, separately, their spectrum (eigenvalues, with the pairs' residuals when computed). KoopmanModel and Eigenfunction objects are built
// from the pair; the expansion itself carries no eigenvalues.
struct EigenfunctionFitResult {
    FunctionExpansion functions;
    Spectrum spectrum;
    EigenfunctionFitDiagnostics diagnostics;
};
// Polynomial eigenfunctions given by their coefficients, without a fit: normalized monomials
// ((q - center) / scale)^exponents (feature x coordinate; empty center/scale mean 0/1), mode-major
// coefficients, on `domain` (default: the physical state of that dimension). Validates the shape,
// finite values, unique monomials, and that no mode is identically zero; no accuracy is invented.
EigenfunctionFitResult polynomialEigenfunctions(std::size_t dimension, std::vector<std::size_t> exponents,
                                                std::vector<std::complex<double>> coefficients, Spectrum spectrum,
                                                std::vector<double> center = {}, std::vector<double> scale = {},
                                                std::optional<FunctionDomain> domain = std::nullopt);
} // namespace dynsyslab::core
