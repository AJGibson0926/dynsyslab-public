#pragma once
// K-ROM archives: the one data-object kind "krom" (core/DataObject.hpp), schema dynsyslab.krom/v1.
//
// A krom is a KoopmanModel: spectrum, eigenfunctions (samples and named representations), and the
// optional input maps B and readouts C. Arrays:
//   spectrum (SpectrumArchive.hpp: eigenvalues, the eigenpairs' residuals, pseudospectrum, continuous
//   part), eigenfunction_values (samples x modes), training_time,
//   training_sample_indices; representations: expansion_coefficients[__<name>] (modes x terms),
//   interpolation_points (samples x domain dimension; its values are eigenfunction_values),
//   window_map (modes x window rows); readouts C: reconstruction_names, reconstruction__<key>
//   (outputs x modes), reconstruction_offset__<key>, reconstruction_error__<key> (rms, max,
//   relative); input maps B: control_<name>_*, discrete_input_map_<i>, input_map_manifest_json;
//   spectral_diagnostics_json (invariance, wave packets, stochastic residuals, RKHS adjoints).
// metadata_json["krom"]: origin, input_channels, representations (the descriptions this archive
// stores), order (every name; the first is the one the model evaluates, also
// "default_representation"), parent.
//
// Parents. A step that adds to a krom it read (B, C, another representation, a mode selection)
// publishes a new krom with parent = {path, sha256, modes}: the parent's modes it keeps, in order. It
// stores only what it adds; spectrum and samples come from the root (a child that attaches residuals
// or a pseudospectrum stores that evidence for its modes). Readers load the parent (refused
// when it is missing or its hash changed), keep those modes, then lay this archive's parts over them
// (a representation, readout or map of the same name replaces the parent's). Archives stay immutable.
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/KoopmanModel.hpp>
#include <dynsyslab/core/ScientificCompatibility.hpp>

#include <filesystem>
#include <map>
#include <string>
#include <vector>

namespace dynsyslab::core {

inline constexpr const char* kKromSchema = "dynsyslab.krom/v1";

// Loads any krom, following its parents. The model evaluates the binding's representation, else the
// archive's first; `selectedModes` keeps a subset. The loaded model's parent is this archive.
[[nodiscard]] KoopmanModel loadKoopmanModel(
    const std::filesystem::path& path, const SystemConfiguration& system,
    const std::vector<std::size_t>& selectedModes = {},
    const ModelBindingOptions& binding = {});
[[nodiscard]] KoopmanModel loadKoopmanModel(const std::filesystem::path& path,
                                            const std::vector<std::size_t>& selectedModes = {});
// The archive's own system, with binding options.
[[nodiscard]] KoopmanModel loadKoopmanModel(const std::filesystem::path& path, const ModelBindingOptions& binding);
// The interpolation, checked to be on the named domain ("state", "input", ...).
[[nodiscard]] KoopmanModel loadKoopmanModelOnDomain(const std::filesystem::path& path,
                                                    const std::string& domain, const json::Value& method = {});
// A krom held in memory (its parents, if any, are read from disk), on the archive's own system. It
// has no parent of its own (it has no path).
[[nodiscard]] KoopmanModel readKoopmanModel(const NpzRecord& archive);

// Whether an archive (or its metadata text) is a krom.
[[nodiscard]] bool isKrom(const NpzRecord& archive);
[[nodiscard]] bool isKromMetadata(const std::string& metadataJson);


// The krom archive of `model`, with `metadataJson` as the base of its metadata (its data_object
// envelope is kept, with kind krom) and `extraArrays` (properties of the producing method). A model
// with a parent stores only its own parts. `inputChannels` describe the record it was fitted on.
[[nodiscard]] NpzRecord makeKoopmanArchive(const KoopmanModel& model,
    const std::string& metadataJson = {}, std::map<std::string, NpyArray> extraArrays = {},
    const std::vector<ChannelDescriptor>& inputChannels = {});
void saveKoopmanModel(const std::filesystem::path& path, const KoopmanModel& model,
                      const std::string& metadataJson = {}, std::map<std::string, NpyArray> extraArrays = {});

// Loads one independently archived constant or bilinear control map. Prefixes support evaluation
// bundles such as "evaluation_0_matrix" without coupling the KROM to the bundle after import.
[[nodiscard]] ApproximateControlMap loadApproximateControlMap(
    const std::filesystem::path& path, const std::string& prefix, std::size_t modes, std::size_t inputs);
// Runtime attachment checks the exact ordered eigenfunctions and time convention. A map that records
// the identity of another K-ROM is refused. A map that records none (one fitted outside the suite)
// is refused unless `externalMap`: then it is attached on the caller's word, and its provenance says
// so (dimension agreement alone is not evidence).
[[nodiscard]] ApproximateControlMap loadApproximateControlMap(
    const std::filesystem::path& path, const std::string& prefix, const KoopmanModel& model,
    bool externalMap = false);

}  // namespace dynsyslab::core
