#pragma once

#include <dynsyslab/core/KoopmanTraining.hpp>
#include <dynsyslab/core/Observable.hpp>

#include <complex>
#include <cstddef>
#include <string>
#include <vector>

namespace dynsyslab::core {

struct EigenvalueOptimizationOptions {
    std::vector<std::vector<std::complex<double>>> initialDiscreteGroups;
    double lowerBound{-1.05};
    double upperBound{1.05};
    std::size_t maximumIterations{500};
    std::size_t memory{8};
    double tolerance{1.0e-10};
    bool preserveConjugatePairs{true};
};

struct EigenvalueOptimizationDiagnostic {
    std::size_t outputIndex{};
    std::size_t modeCount{};
    std::size_t iterations{};
    std::size_t functionEvaluations{};
    double normalizedSquaredResidual{};
    double relativeResidual{};
    bool converged{};
    std::string message;
};

struct EigenvalueOptimizationResult {
    std::vector<std::complex<double>> discreteEigenvalues;
    std::vector<std::complex<double>> continuousEigenvalues;
    std::vector<std::vector<std::size_t>> modeGroups;
    std::vector<EigenvalueOptimizationDiagnostic> diagnostics;
    double deltaTime{};
};

// Variable-projection eigenvalue fit: trajectory boundary amplitudes are
// eliminated analytically at every bounded L-BFGS objective evaluation. This
// stage returns a spectrum only and never constructs eigenfunctions.
[[nodiscard]] EigenvalueOptimizationResult optimizeTrajectoryEigenvalues(
    const KoopmanTrajectoryEnsemble& ensemble,
    const Observable& observable,
    const EigenvalueOptimizationOptions& options);

}  // namespace dynsyslab::core
