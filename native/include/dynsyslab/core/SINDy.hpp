#pragma once

#include <complex>
#include <cstddef>
#include <string>
#include <vector>

#include <dynsyslab/core/Observable.hpp>

namespace dynsyslab::core {

struct KoopmanSINDyOptions {
    std::size_t polynomialDegree{4};
    double threshold{0.05};
    std::size_t maximumIterations{20};
    bool normalizeStates{true};
    bool includeConstant{true};
    bool reuseConjugateFits{true};
    bool materializeTrainingValues{false};
};

// Compatibility fit DTO. New runtime code uses PolynomialEigenfunctionSet;
// new SINDy callers use EigenfunctionFitResult in EigenfunctionFit.hpp.
struct KoopmanSINDyResult {
    // Construct the SAME polynomial representation without fitting. Coefficients
    // are mode-major; exponents are feature-major (arbitrary unique monomials).
    // Eigenvalues are continuous-time rates, never sampled multipliers. Empty
    // center/scale mean raw variables (0/1). No accuracy metrics are fabricated.
    [[nodiscard]] static KoopmanSINDyResult fromCoefficients(
        std::size_t coordinateDimension,
        std::vector<std::size_t> exponents,
        std::vector<std::complex<double>> coefficients,
        std::vector<std::complex<double>> continuousEigenvalues,
        std::vector<double> center = {}, std::vector<double> scale = {});

    [[nodiscard]] static KoopmanSINDyResult fromDiscreteCoefficients(
        std::size_t coordinateDimension, std::vector<std::size_t> exponents,
        std::vector<std::complex<double>> coefficients,
        std::vector<std::complex<double>> discreteEigenvalues, double deltaTime,
        std::vector<double> center = {}, std::vector<double> scale = {},
        std::vector<int> logarithmBranches = {});
    // dt=0 uses the saved sampling interval; a missing interval is an error.
    // Changing dt changes the multipliers, not the stored continuous rates.
    [[nodiscard]] std::vector<std::complex<double>> discreteEigenvalues(double dt = 0) const;

    std::vector<std::complex<double>> continuousEigenvalues;
    std::vector<std::complex<double>> coefficients; // eigenvalue x feature
    std::vector<std::complex<double>> trainingValues; // sample x eigenvalue, optional
    std::vector<std::size_t> exponents; // feature x state
    std::vector<double> center;
    std::vector<double> scale;
    std::vector<double> residualRms;
    std::vector<double> residualMaximum;
    std::size_t samples{};
    std::size_t stateDimension{};
    std::size_t dynamicalStateDimension{};
    std::size_t featureCount{};
    std::size_t directFits{};
    std::vector<std::string> coordinateObservableNames;
    std::string coordinateObservableDefinition; // versioned exact expression definition, when available
    // Empty or "instantaneous-observable" means y=h(x).  "causal-history"
    // means the polynomial variables are supplied by a separate stateful
    // history projector and cannot be inferred from the instantaneous state.
    std::string coordinateProvider;
    double sampleTime{}; // output sampling interval, NOT an adaptive integrator step
    std::string eigenvalueOrigin{"continuous-time"};
    std::vector<int> logarithmBranches; // per mode for discrete-origin rates
};

// Canonical graded monomial order shared by fitting and archive import.
[[nodiscard]] std::vector<std::size_t> polynomialExponentTable(
    std::size_t dimension, std::size_t degree, bool includeConstant);

// Fits polynomial generator-PDE eigenfunctions for an explicitly supplied
// continuous-time spectrum. Velocities may be supplied from the exact ODE;
// otherwise second-order finite differences are used trajectory-locally.
[[nodiscard]] KoopmanSINDyResult fitKoopmanEigenfunctions(
    const std::vector<double>& states, std::size_t samples,
    std::size_t stateDimension, const std::vector<double>& time,
    const std::vector<std::complex<double>>& continuousEigenvalues,
    const std::vector<double>& velocities = {},
    const KoopmanSINDyOptions& options = {});

// Fits eigenfunctions of y=h(x). Exact state velocities are transformed by
// ydot=J_h(x)xdot before constructing the generator-PDE library.
[[nodiscard]] KoopmanSINDyResult fitKoopmanEigenfunctionsInCoordinates(
    const std::vector<double>& dynamicalStates, std::size_t samples,
    std::size_t dynamicalStateDimension, const std::vector<double>& time,
    const std::vector<std::complex<double>>& continuousEigenvalues,
    const Observable& coordinates,
    const std::vector<double>& dynamicalVelocities = {},
    const KoopmanSINDyOptions& options = {});

[[nodiscard]] std::complex<double> evaluatePolynomialEigenfunction(
    const KoopmanSINDyResult& model, std::size_t eigenfunction,
    const double* state);

[[nodiscard]] std::vector<std::complex<double>> gradientPolynomialEigenfunction(
    const KoopmanSINDyResult& model, std::size_t eigenfunction,
    const double* state);

}  // namespace dynsyslab::core
