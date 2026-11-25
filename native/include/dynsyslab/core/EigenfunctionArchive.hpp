#pragma once
// Archive contract for EigenfunctionSet: dynsyslab.eigenfunction-objects/v1.
//
// Arrays: continuous_eigenvalues, discrete_eigenvalues, eigenfunction_values (samples x modes),
// training_time, training_sample_indices (rows of the input series the origin points to),
// interpolation_points (the interpolation's base points q_k, samples x domain dimension),
// window_map (modes x rows: row j is eigenfunction j's own linear map from a raw window).
// metadata_json["eigenfunction_objects"]: origin, input channels, representations (in order); the
// interpolation's description holds its domain (with its state map) and its interpolant method.
// No input values are copied: the window map reads them from the input series by its layout.
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/HankelRepresentations.hpp>

namespace dynsyslab::core {

inline constexpr const char* kEigenfunctionObjectsSchema = "dynsyslab.eigenfunction-objects/v1";

// Adds the set's arrays to `archive` and returns the metadata block for "eigenfunction_objects".
// `modes` (Hankel modes, rows x rank) are written by the caller as the existing "modes" array.
[[nodiscard]] json::Value writeEigenfunctionSet(NpzRecord& archive, const EigenfunctionSet& set,
                                                const std::vector<ChannelDescriptor>& inputChannels);

[[nodiscard]] bool hasEigenfunctionObjects(const NpzRecord& archive);

// Rebuilds the set. `kinds`, when given, loads only those representations (skipping e.g. a large
// window map).
[[nodiscard]] EigenfunctionSet readEigenfunctionSet(const NpzRecord& archive, const std::vector<std::string>& kinds = {});

// One eigenfunction on its own: its eigenvalue and the chosen representation (default: the first),
// independent of the method that produced it. Its window map reads raw windows of the channels and
// layout its description names (HankelWindowLayout::window_at), so new windows of the same length
// evaluate without the producing method.
[[nodiscard]] Eigenfunction readEigenfunction(const NpzRecord& archive, std::size_t mode,
                                             const std::string& representation = {},
                                             const ObservableCatalog* catalog = nullptr);


}  // namespace dynsyslab::core
