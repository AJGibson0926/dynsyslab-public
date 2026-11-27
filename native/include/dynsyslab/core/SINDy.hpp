#pragma once

#include <complex>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

#include <dynsyslab/core/EigenfunctionFit.hpp>
#include <dynsyslab/core/ExpansionDictionary.hpp>
#include <dynsyslab/core/Observable.hpp>

namespace dynsyslab::core {

struct KoopmanSINDyOptions {
    std::size_t polynomialDegree{4};
    double threshold{0.05};
    std::size_t maximumIterations{20};
    bool normalizeStates{true};
    bool includeConstant{true};
    bool reuseConjugateFits{true};
    // Solve each thresholded null vector from feature-sized Gram blocks formed
    // once (ThetaᵀTheta, GammaᵀTheta, GammaᵀGamma) instead of a tall SVD per
    // eigenvalue and iteration. A mode whose Gram system is too ill-conditioned
    // for its smallest eigenvector to be trusted falls back to the tall SVD.
    bool gramNullspace{true};
    // Anchored fit (optional). With reference values psi_j(k) of every eigenfunction at anchor states
    // x_k, each coefficient vector solves
    //   min_xi ||Theta(X_a) xi - psi_j||^2 + w ||(Gamma - lambda_j Theta) xi||^2
    // instead of the null space of the generator equation alone. The reference fixes the free scale
    // and phase and keeps the fit on the eigenfunction the values describe. `anchorWeight` is relative:
    // w = anchorWeight * trace(Theta_aᵀTheta_a) / trace(G(lambda_j)). States are in the same variables
    // as the fit's (dynamical states for fitKoopmanEigenfunctionsInCoordinates).
    const std::vector<double>* anchorStates{nullptr};                 // anchors x dimension
    const std::vector<std::complex<double>>* anchorValues{nullptr};   // anchors x eigenvalues
    double anchorWeight{1.0};
    // Positive per-component multipliers after normalization, q' = (q - center) / scale * m (empty:
    // none; one value applies to every component).
    std::vector<double> rescale;
    // Sample counts of the trajectories the rows pool (empty: one trajectory). Time must increase and
    // finite differences are taken within each trajectory, so no difference crosses a boundary. Rows
    // with supplied velocities and no trajectory (a state-space domain) are segments of one sample.
    std::vector<std::size_t> segments;
    // How each null vector is sparsified: "iterative-relative" (sequential thresholding of the
    // column-normalized generator, relative to its largest entry; the default) or "single-absolute"
    // (the smallest right singular vector of the unnormalized generator once, then every entry below
    // `threshold` in magnitude zeroed, with no refit).
    std::string sparsification{"iterative-relative"};
    // Gamma = grad Theta(q) . qdot from the dictionary gradient ("analytic", exact for monomials) or
    // by the directional central difference (Theta(q + h qdot) - Theta(q - h qdot)) / 2h with
    // h = directionalStep * max(1, max_i |q_i|) ("finite-difference").
    std::string directionalDerivative{"analytic"};
    double directionalStep{1e-6};
    // Fit in this dictionary instead of the polynomial one: any ExpansionDictionary terms (radial,
    // trigonometric, expression, ...), each evaluated on the raw coordinates with its own
    // normalization, so normalizeStates and rescale do not apply.
    std::shared_ptr<const ExpansionDictionary> dictionary;
};

// Canonical graded monomial order shared by fitting and archive import.
[[nodiscard]] std::vector<std::size_t> polynomialExponentTable(
    std::size_t dimension, std::size_t degree, bool includeConstant);

// Fits polynomial generator-PDE eigenfunctions for an explicitly supplied
// continuous-time spectrum. Velocities may be supplied from the exact ODE;
// otherwise second-order finite differences are used trajectory-locally.
// The result is on the physical state of that dimension (fitSindyEigenfunctions names
// another domain); its diagnostics are the generator residuals Dpsi(x) xdot - lambda psi(x).
[[nodiscard]] EigenfunctionFitResult fitKoopmanEigenfunctions(
    const std::vector<double>& states, std::size_t samples,
    std::size_t stateDimension, const std::vector<double>& time,
    const std::vector<std::complex<double>>& continuousEigenvalues,
    const std::vector<double>& velocities = {},
    const KoopmanSINDyOptions& options = {});

// Fits eigenfunctions of y=h(x). Exact state velocities are transformed by
// ydot=J_h(x)xdot before constructing the generator-PDE library.
[[nodiscard]] EigenfunctionFitResult fitKoopmanEigenfunctionsInCoordinates(
    const std::vector<double>& dynamicalStates, std::size_t samples,
    std::size_t dynamicalStateDimension, const std::vector<double>& time,
    const std::vector<std::complex<double>>& continuousEigenvalues,
    const Observable& coordinates,
    const std::vector<double>& dynamicalVelocities = {},
    const KoopmanSINDyOptions& options = {});

// Koopman eigenvalues from data alone, in the same polynomial library as the fits.
// Kaiser, Kutz & Brunton, Mach. Learn.: Sci. Technol. 2:035023 (2021), Sec. 5.1, Algorithm 1:
//   1. initial guesses: the eigenvalues of K = Theta^+ Gamma, the generator in the library
//      (generator EDMD: Klus et al., Physica D 406:132416, 2020), formed from the Gram blocks as
//      K = (ThetaᵀTheta)^-1 ThetaᵀGamma;
//   2. for each guess, the sparse null vector xi of (Gamma - lambda Theta) (here by sequential
//      thresholding, as in the fits, in place of their alternating-directions sparse-vector search);
//   3. the Rayleigh update lambda <- xi^H K xi / xi^H xi, repeated until lambda settles;
//   4. validation along the training series (their Eq. 47): the relative error of
//      phi(x(t)) against c e^{lambda t}, with c fitted by least squares rather than taken from
//      phi(x(0)). Eigenvalues above `validationTolerance` and the trivial lambda = 0 are dropped.
// The training series must be one trajectory; a finite library may not be closed under the dynamics,
// so the estimate can miss eigenvalues or contain spurious ones (the validation catches the latter).
struct KoopmanSINDyEigenvalueOptions {
    std::size_t refinementIterations{50};
    double refinementTolerance{1e-10};  // relative change of lambda
    double validationTolerance{0.1};
    bool keepTrivial{false};
};
struct KoopmanSINDyEigenvalueEstimate {
    std::vector<std::complex<double>> eigenvalues;  // continuous-time, kept
    std::vector<std::complex<double>> initial;      // their library-generator guesses
    std::vector<double> validationError;            // per kept eigenvalue
    std::vector<std::size_t> iterations;            // refinement steps per kept eigenvalue
    std::size_t candidates{}, rejected{};
};
[[nodiscard]] KoopmanSINDyEigenvalueEstimate estimateKoopmanEigenvaluesInCoordinates(
    const std::vector<double>& dynamicalStates, std::size_t samples,
    std::size_t dynamicalStateDimension, const std::vector<double>& time,
    const Observable& coordinates,
    const std::vector<double>& dynamicalVelocities = {},
    const KoopmanSINDyOptions& options = {},
    const KoopmanSINDyEigenvalueOptions& estimate = {});

// The same generator fit on explicit coordinates `values` (samples x coordinates.dimension), with
// diagnostics kept separate. Supplied derivatives may be measured or modeled.
[[nodiscard]] EigenfunctionFitResult fitSindyEigenfunctions(
    const std::vector<double>& values, std::size_t samples, const std::vector<double>& time,
    const std::vector<std::complex<double>>& continuousRates, FunctionDomain coordinates,
    const std::vector<double>& derivatives = {}, const KoopmanSINDyOptions& options = {});

}  // namespace dynsyslab::core
