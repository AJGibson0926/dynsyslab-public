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

}  // namespace dynsyslab::core::detail
