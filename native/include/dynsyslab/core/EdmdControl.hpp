#pragma once

#include <dynsyslab/core/DmdControl.hpp>
#include <dynsyslab/core/SINDy.hpp>

namespace dynsyslab::core {

struct EdmdControlOptions {
    DmdControlOptions regression;
    std::size_t polynomialDegree{4};
    bool includeConstant{true};
    bool normalize{true};
};

// Same total-degree monomial order and population-standard-deviation scaling
// as SINDy. Normalization is fitted on X only, then frozen for Y and deployment.
struct EdmdControlDictionary {
    std::size_t dimension{}, degree{4};
    bool includeConstant{true};
    std::vector<std::size_t> exponents; // feature x supplied coordinate
    std::vector<double> center, scale;
};

struct EdmdControlResult {
    EdmdControlDictionary dictionary;
    DmdControlResult lifted; // q here is Theta(original snapshot), not raw state.
    double snapshotOneStepRelativeError{};
    // Explicit export to the existing polynomial eigenfunction class, using
    // principal log(mu)/dt rates. Throws for zero multipliers or an unusable
    // eigenbasis. This does NOT convert the discrete input map to continuous B.
    [[nodiscard]] EigenfunctionFitResult eigenfunctionSet() const;
    [[nodiscard]] KoopmanSINDyResult expansion() const;
};

// Pure measured-data API. No derivatives, EOM, history construction, or sparse
// SINDy optimization. Inputs remain unlifted; this fits an additive discrete B.
// Explicit sample-major pairs can be pooled without joining trajectory endpoints.
[[nodiscard]] EdmdControlResult extendedDynamicModeDecompositionWithControl(
    const std::vector<double>& x, const std::vector<double>& y,
    const std::vector<double>& inputs, std::size_t pairs,
    std::size_t dimension, std::size_t inputDimension, double deltaTime,
    const EdmdControlOptions& options = {});

[[nodiscard]] std::vector<double> liftEdmdControlSnapshots(
    const EdmdControlDictionary& dictionary, const std::vector<double>& snapshots,
    std::size_t samples);
// Polynomial coefficients are lifted.eigenfunctionCoefficients; evaluate them
// against the saved normalized monomials, NEVER against raw coordinates.
[[nodiscard]] std::vector<std::complex<double>> evaluateEdmdControlEigenfunctions(
    const EdmdControlResult& model, const std::vector<double>& snapshots,
    std::size_t samples);
// Returns predicted snapshots in original stored units. Uses the degree-one
// dictionary entries and fixed affine inverse normalization (not a new fit).
[[nodiscard]] std::vector<double> predictEdmdControl(
    const EdmdControlResult& model, const std::vector<double>& snapshots,
    const std::vector<double>& inputs, std::size_t samples);

} // namespace dynsyslab::core
