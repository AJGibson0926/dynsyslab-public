#pragma once

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/KoopmanModel.hpp>
#include <dynsyslab/core/ScientificCompatibility.hpp>

#include <filesystem>
#include <string>
#include <vector>

namespace dynsyslab::core {
// Scattered eigenfunction objects bound on a named interpolation domain ("state", "input", ...).
[[nodiscard]] KoopmanModel loadKoopmanModelOnDomain(const std::filesystem::path& path,
                                                    const std::string& domain, const json::Value& method = {});


// Loads either an analytic-polynomial or scattered-interpolation KROM. This is
// the normal runtime boundary for controllers, validation, and Viewer.
[[nodiscard]] KoopmanModel loadKoopmanModel(
    const std::filesystem::path& path, const SystemConfiguration& system,
    const std::vector<std::size_t>& selectedModes = {},
    const ModelBindingOptions& binding = {});
[[nodiscard]] KoopmanModel loadKoopmanModel(
    const std::filesystem::path& path,
    const std::vector<std::size_t>& selectedModes = {});

// The archive's own system, with binding options (e.g. reconstructions objects to attach).
[[nodiscard]] KoopmanModel loadKoopmanModel(const std::filesystem::path& path, const ModelBindingOptions& binding);

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

// A reconstructions object (core/DataObject.hpp kind "reconstructions"): observables written over one
// eigenfunctions object, g ~ C psi (+ offset). Arrays: reconstruction_names (one per line),
// reconstruction__<key> (outputs x modes), optional reconstruction_offset__<key> and
// reconstruction_error__<key> (rms, max, relative). metadata_json["reconstructions"] names the
// basis: {"basis": {"path", "sha256"}, "modes"}.
[[nodiscard]] json::Value reconstructionsBlock(const std::filesystem::path& eigenfunctions, std::size_t modes);
// Adds the reconstructions' arrays to an in-memory archive of the eigenfunctions they are written over,
// after checking that `eigenfunctionsPath` is their basis (by SHA-256) and the mode count agrees.
void mergeReconstructions(NpzRecord& eigenfunctions, const std::filesystem::path& eigenfunctionsPath,
                          const NpzRecord& reconstructions);

void saveKoopmanModel(const std::filesystem::path& path,
                      const KoopmanModel& model,
                      const std::string& metadataJson = {},
                      std::map<std::string, NpyArray> extraArrays = {});

}  // namespace dynsyslab::core
