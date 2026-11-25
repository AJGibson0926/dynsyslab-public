#pragma once

#include <dynsyslab/core/Observable.hpp>
#include <dynsyslab/core/KoopmanTraining.hpp>
#include <dynsyslab/core/Reconstruction.hpp>
#include <dynsyslab/core/FunctionInterpolation.hpp>

#include <complex>
#include <cstddef>
#include <vector>

namespace dynsyslab::core {

struct BoundaryEigenfunctionOptions {
    std::vector<std::vector<std::size_t>> modeGroups;
    std::size_t interpolationNeighbors{8};
    double gradientRelativeStep{1.0e-6};
    double singularValueTolerance{1.0e-12};
};

struct BoundaryEigenfunctionResult {
    InterpolationData eigenfunctions;                        // the interpolation representation
    std::vector<std::complex<double>> continuousEigenvalues;  // their eigenvalues (the supplied ones)
    ReconstructionResult observableReconstruction;
    std::vector<double> relativeResidualByOutput;
    std::vector<std::vector<std::size_t>> modeGroups;
};

// Regresses trajectory-specific boundary amplitudes for immutable supplied
// eigenvalues, then materializes their exact temporal powers as scattered
// state-space eigenfunctions. Eigenvalue extraction is deliberately separate.
[[nodiscard]] BoundaryEigenfunctionResult fitBoundaryEigenfunctions(
    const std::vector<std::complex<double>>& continuousEigenvalues,
    double deltaTime,
    const KoopmanTrajectoryEnsemble& ensemble,
    const Observable& observable,
    const BoundaryEigenfunctionOptions& options = {});

}  // namespace dynsyslab::core
