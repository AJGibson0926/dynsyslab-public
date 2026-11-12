#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

#include <dynsyslab/core/Spectral.hpp>

namespace dynsyslab::core {

struct PCAOptions {
    std::size_t rank{}; // zero retains every singular direction
    bool center{true};
};

struct PCAResult {
    std::vector<double> mean;
    std::vector<double> components; // rank x feature
    std::vector<double> singularValues;
    std::vector<double> scores; // sample x rank
    std::vector<double> explainedVarianceRatio;
    std::size_t samples{};
    std::size_t features{};
    std::size_t rank{};
};

[[nodiscard]] PCAResult principalComponentAnalysis(
    const std::vector<double>& values, std::size_t samples,
    std::size_t features, const PCAOptions& options = {});

struct BaggedPCAOptions : PCAOptions {
    std::size_t bags{32};
    std::optional<std::uint64_t> seed;
};

struct BaggedPCAResult {
    PCAResult decomposition;
    std::size_t bags{};
    std::uint64_t seed{};
};

[[nodiscard]] BaggedPCAResult baggedPrincipalComponentAnalysis(
    const std::vector<double>& values, std::size_t samples,
    std::size_t features, const BaggedPCAOptions& options = {});

struct RobustPCAOptions {
    double lambda{}; // zero selects 1/sqrt(max(samples, features))
    double tolerance{1.0e-7};
    std::size_t maximumIterations{1000};
};

struct RobustPCAResult {
    std::vector<double> lowRank;
    std::vector<double> sparse;
    std::size_t samples{};
    std::size_t features{};
    std::size_t iterations{};
    double relativeResidual{};
    double lambda{};
    bool converged{};
};

[[nodiscard]] RobustPCAResult robustPrincipalComponentAnalysis(
    const std::vector<double>& values, std::size_t samples,
    std::size_t features, const RobustPCAOptions& options = {});

struct HankelPODOptions {
    std::size_t matrixSize{10};
    std::size_t delayStride{1};
    std::size_t rank{};
    bool center{};
};

struct HankelPODResult {
    SVDResult decomposition;
    std::vector<double> mean;
    std::size_t embeddedSamples{};
    std::size_t embeddedFeatures{};
    std::size_t matrixSize{};
    std::size_t delayStride{};
};

[[nodiscard]] HankelPODResult fitHankelPOD(const HankelMatrix& matrix, std::size_t rank, bool center);

[[nodiscard]] HankelPODResult hankelProperOrthogonalDecomposition(
    const std::vector<double>& values, std::size_t samples,
    std::size_t features, const HankelPODOptions& options = {});

}  // namespace dynsyslab::core
