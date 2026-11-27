#pragma once

// Output-grid array operations for ODE solution archives: truncation to computed rows,
// post-solve resampling with raw_* copies, and checkpoint (partial) archives.

#include "SolveSettings.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/workflow/Json.hpp>

#include <cstddef>
#include <filesystem>
#include <map>
#include <string>
#include <vector>

namespace dynsyslab::workflow::detail {

// Keep the first `rows` rows of every array whose leading dimension is `sampleCount`.
void truncateOutputRows(core::NpzRecord& archive, std::size_t sampleCount, std::size_t rows);

// Resample every real output-grid array (leading dimension = time.size()) onto
// settings.times; the source keeps only its first `rows` rows (a truncated trajectory gives
// NaN after its last computed time). The original arrays are stored as raw_<name> and
// raw_time. Returns the metadata block recorded as `interpolation`.
json::Value resampleOutputArrays(core::NpzRecord& archive, const std::vector<double>& time,
                                 std::size_t rows, const ResampleSettings& settings);

// Partial archive written at a checkpoint: `<archive>.partial.npz` beside the final archive.
[[nodiscard]] std::filesystem::path partialArchivePath(const std::filesystem::path& archive);
void writePartialArchive(const std::filesystem::path& path, const core::IntegrationCheckpoint& point,
                         const std::string& metadataJson);

// Partial archive of a delay or lattice run: time[0..rows), values rows x width, optional
// int64 arrays (for example the current lattice), metadata. Written atomically.
void writePartialSeries(const std::filesystem::path& path, const double* time, std::size_t rows,
                        const double* values, std::size_t width, const std::string& metadataJson,
                        const std::map<std::string, core::NpyArray>& extra = {});

// history.checkpoint_interval of a solve job in its own clock units (0: no checkpoints).
[[nodiscard]] double checkpointInterval(const json::Value& job);

}  // namespace dynsyslab::workflow::detail
