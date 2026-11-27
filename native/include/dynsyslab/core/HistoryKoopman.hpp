#pragma once

#include <dynsyslab/core/EigenfunctionFit.hpp>
#include <dynsyslab/core/SINDy.hpp>
#include <dynsyslab/core/HistoryObservation.hpp>

#include <complex>
#include <cstddef>
#include <optional>
#include <vector>

namespace dynsyslab::core {

// A fixed scalar-observation delay basis fitted across independent
// trajectories. No Hankel window is ever allowed to cross a trajectory
// boundary.
struct PooledHankelOptions {
    std::size_t delays{100};
    std::size_t rank{4};
    std::size_t oversampling{4};
    std::size_t subspaceIterations{8};
    std::size_t seed{20261002};
    bool normalizeByRootMeanSquare{true};
};


[[nodiscard]] HankelLeftBasis fitPooledHankelBasis(
    const std::vector<std::vector<double>>& observations,
    const PooledHankelOptions& options = {});


struct PooledPolynomialEDMDOptions {
    std::size_t degree{5};
    std::size_t pairStride{1};
    double sampleTime{1.0};
    double relativeRidge{1.0e-10};
    bool includeConstant{true};
};

struct PooledPolynomialEDMD {
    PolynomialDictionary basis;  // normalized monomials of the history coordinates xi1, xi2, ...
    std::vector<std::complex<double>> coefficients;  // mode x feature, one row per discrete eigenvalue
    std::size_t samples{};
    std::vector<std::complex<double>> discreteEigenvalues;
    std::vector<double> operatorMatrix; // feature x feature, row-major
    std::vector<double> featureScale;
    std::size_t trajectories{};
    std::size_t pairs{};
    std::size_t pairStride{};
    double pairTime{};
    double ridge{};
};

// Fits one discrete polynomial Koopman operator from independent coordinate
// runs. Snapshot pairs are formed trajectory-locally and never bridge runs.
[[nodiscard]] PooledPolynomialEDMD fitPooledPolynomialEDMD(
    const std::vector<std::vector<double>>& coordinateRuns,
    std::size_t coordinateDimension,
    const PooledPolynomialEDMDOptions& options = {});

struct EDMDModeSearchOptions {
    double minimumDiscreteMagnitude{0.90};
    double maximumDiscreteMagnitude{1.005};
    double minimumContinuousReal{-1.0};
    double maximumContinuousReal{0.1};
    double minimumContinuousImaginary{4.0};
    double maximumContinuousImaginary{7.0};
    bool positiveImaginaryOnly{true};
};

struct EDMDModeCandidate {
    std::size_t index{};
    std::complex<double> discreteEigenvalue;
    std::complex<double> continuousEigenvalue;
    std::vector<std::complex<double>> coefficients; // canonical raw monomials
    double meanPairResidual{};
    double minimumPhaseRSquared{};
};

[[nodiscard]] std::vector<EDMDModeCandidate> selectPolynomialEDMDModes(
    const PooledPolynomialEDMD& model,
    const std::vector<std::vector<double>>& coordinateRuns,
    const EDMDModeSearchOptions& options = {});

// Extracts one selected mode into the canonical polynomial eigenfunction
// contract used by KoopmanModel, archive I/O, validation, reconstruction and
// control. The history transform remains a separate, stateful observation
// provider.
[[nodiscard]] EigenfunctionFitResult selectedPolynomialEigenfunction(
    const PooledPolynomialEDMD& model,
    const EDMDModeCandidate& candidate);

}  // namespace dynsyslab::core
