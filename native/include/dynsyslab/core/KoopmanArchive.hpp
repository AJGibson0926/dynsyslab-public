#pragma once

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/KoopmanModel.hpp>
#include <dynsyslab/core/ScientificCompatibility.hpp>

#include <filesystem>
#include <string>
#include <vector>

namespace dynsyslab::core {

// Loads either an analytic-polynomial or scattered-interpolation KROM. This is
// the normal runtime boundary for controllers, validation, and Viewer.
[[nodiscard]] KoopmanModel loadKoopmanModel(
    const std::filesystem::path& path, const SystemConfiguration& system,
    const std::vector<std::size_t>& selectedModes = {},
    const ModelBindingOptions& binding = {});
[[nodiscard]] KoopmanModel loadKoopmanModel(
    const std::filesystem::path& path,
    const std::vector<std::size_t>& selectedModes = {});

// Evaluator-only import preserving every explicitly stored legacy conjugate mode.
// Control loaders continue to retain the independent split-real realization.
[[nodiscard]] KoopmanModel loadKoopmanEigenfunctionModel(
    const std::filesystem::path& path, const SystemConfiguration& system,
    const ModelBindingOptions& binding = {});
[[nodiscard]] KoopmanModel loadKoopmanEigenfunctionModel(const std::filesystem::path& path);

// Direct native archive boundary. The system configuration is explicit so a
// KROM cannot silently bind itself to default physics.
[[nodiscard]] KoopmanModel loadPolynomialKoopmanModel(
    const std::filesystem::path& path, const SystemConfiguration& system,
    const std::vector<std::size_t>& selectedModes = {},
    const ModelBindingOptions& binding = {});

// Loads a native self-contained KROM using the resolved system configuration
// embedded by savePolynomialKoopmanModel.  Legacy archives still require the
// explicit-system overload above.
[[nodiscard]] KoopmanModel loadPolynomialKoopmanModel(
    const std::filesystem::path& path,
    const std::vector<std::size_t>& selectedModes = {});

// Loads one independently archived constant or bilinear control map. Prefixes
// support evaluation bundles such as "evaluation_0_matrix" without coupling
// the KROM to the bundle after import.
[[nodiscard]] ApproximateControlMap loadApproximateControlMap(
    const std::filesystem::path& path, const std::string& prefix,
    std::size_t modes, std::size_t inputs);
// Runtime attachment checks the exact ordered eigenfunctions and time convention.
// Old unbound maps need an explicit reason; dimension agreement alone is not evidence.
[[nodiscard]] ApproximateControlMap loadApproximateControlMap(
    const std::filesystem::path& path, const std::string& prefix,
    const KoopmanModel& model, const std::string& unverifiedLegacyReason = {});

void savePolynomialKoopmanModel(const std::filesystem::path& path,
                                const KoopmanModel& model,
                                const std::string& metadataJson = {},
                                std::map<std::string, NpyArray> extraArrays = {});

// Canonical encoding for workflow publishers that allocate immutable destinations.
[[nodiscard]] NpzRecord makeKoopmanArchive(const KoopmanModel& model,
    const std::string& metadataJson = {}, std::map<std::string,NpyArray> extraArrays = {});

void saveKoopmanModel(const std::filesystem::path& path,
                      const KoopmanModel& model,
                      const std::string& metadataJson = {},
                      std::map<std::string, NpyArray> extraArrays = {});

}  // namespace dynsyslab::core
