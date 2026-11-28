#pragma once
// The parts of a krom archive (KoopmanArchive.hpp): the eigenfunctions' samples and representations,
// and the readouts C. Each writer stores only what the model owns (KoopmanModel::owns*).
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/KoopmanModel.hpp>

#include <memory>
#include <string>
#include <vector>

namespace dynsyslab::core::detail {

// eigenfunction_values (samples x modes), training_time, training_sample_indices.
void writeKromSamples(NpzRecord& archive, const EigenfunctionSamples& samples);
[[nodiscard]] EigenfunctionSamples readKromSamples(const NpzRecord& archive, std::size_t modes);

// The arrays and descriptions of the representations the model owns (every one when `ownOnly` is false).
[[nodiscard]] json::Value::Array writeKromRepresentations(NpzRecord& archive, const KoopmanModel& model, bool ownOnly);
// Rebuilds stored representations; an interpolation takes its values from `samples`.
[[nodiscard]] std::vector<std::shared_ptr<const FunctionRepresentation>> readKromRepresentations(
    const json::Value& descriptions, const NpzRecord& archive, const EigenfunctionSamples& samples, std::size_t modes);

// reconstruction_names, reconstruction__<key>, reconstruction_offset__<key>, reconstruction_error__<key>.
void writeKromReadouts(NpzRecord& archive, const KoopmanModel& model, bool ownOnly);
void readKromReadouts(const NpzRecord& archive, KoopmanModel& model);

// The compressed operator (compression_basis, compression_action, compression_readout__<key>) and its
// description for the krom block (KromSpectralParts.cpp).
[[nodiscard]] json::Value writeKromCompression(NpzRecord& archive, const KoopmanCompression& compression);
[[nodiscard]] KoopmanCompression readKromCompression(const json::Value& block, const NpzRecord& archive);
// The functions on spectral sets the model owns (spectral_function__<key>) and their descriptions.
[[nodiscard]] json::Value::Array writeKromSpectralFunctions(NpzRecord& archive, const KoopmanModel& model, bool ownOnly);
void readKromSpectralFunctions(const json::Value& descriptions, const NpzRecord& archive, KoopmanModel& model);

}  // namespace dynsyslab::core::detail
