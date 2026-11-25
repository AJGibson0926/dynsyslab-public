#pragma once
// Projection of observables onto eigenfunctions through their shared dictionary, and the closed
// form of an expansion as an expression Observable. (An expansion becomes an ordinary Observable
// through Observable::fromRepresentation; this header adds what only expansions can do.)
#include <dynsyslab/core/Eigenfunction.hpp>
#include <dynsyslab/core/FunctionExpansion.hpp>

#include <optional>

namespace dynsyslab::core {

// The outputs as closed-form real expressions over the domain coordinates (real, imaginary part per
// output), when every term has one. Gives partial, flow and input-flow derivatives and an
// expression definition, exactly as for any expression Observable.
[[nodiscard]] std::optional<Observable> expressionObservable(const FunctionExpansion& expansion,
                                                             ObservableDescriptor descriptor);

// Several expansions on one domain as one multi-output expansion in the union of their dictionaries.
[[nodiscard]] FunctionExpansion stackExpansions(const std::vector<FunctionExpansion>& expansions);

// g_r ~ sum_j B_rj psi_j. "coefficients": least squares on the canonical dictionary coefficients;
// exact (an identity of functions) when every g_r lies in span{psi_j}, giving analytic B, e.g. the
// K-ROM reconstruction of g. "samples": least squares over the supplied points (rows of q).
struct ExpansionProjection {
    std::string innerProduct;
    std::vector<std::complex<double>> coefficients;  // outputs x modes (B)
    std::vector<double> relativeResidual;            // per output, in the chosen inner product
    bool exact{};
    FunctionExpansion residual;                      // g - B psi, canonical dictionary
};
[[nodiscard]] ExpansionProjection projectExpansion(const FunctionExpansion& observable, const FunctionExpansion& onto,
                                                   const std::vector<double>& samplePoints = {}, double tolerance = 1e-10);
// The same for Observable and Eigenfunction objects whose representation is an expansion.
[[nodiscard]] ExpansionProjection projectObservable(const Observable& observable, const std::vector<Eigenfunction>& onto,
                                                    const std::vector<double>& samplePoints = {}, double tolerance = 1e-10);

}  // namespace dynsyslab::core
