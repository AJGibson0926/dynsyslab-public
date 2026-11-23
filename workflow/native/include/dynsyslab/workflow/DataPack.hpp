#pragma once

#include <dynsyslab/workflow/Json.hpp>

#include <filesystem>
#include <string>

namespace dynsyslab::workflow {

// A data pack is a folder of ordinary job files plus a manifest that pins them to one release:
//   pack.json (dynsyslab.data-pack/v1), pack.sha256, README.md, jobs/00-reproduce-all.job.json, ...
// Execution belongs to the job system (kind "job_group"); this layer checks pins and files,
// compares a finished run with the sealed digests, and seals and exports packs.

inline constexpr const char* kDataPackSchema = "dynsyslab.data-pack/v1";

// For a job naming "data_pack": refuses a release other than the pinned one, unless
// runtime.allow_version_mismatch is set, in which case the run is labelled non-reproducing.
// Returns the pin record stored in the ledger (null for a job without a pack).
[[nodiscard]] json::Value checkDataPackPin(const json::Value& job);

// Reads pack.json from a folder. Any dynsyslab.data-pack/vN stays inspectable; only versions
// this release understands may run.
[[nodiscard]] json::Value readDataPack(const std::filesystem::path& folder);
[[nodiscard]] json::Value inspectDataPack(const std::filesystem::path& folder);

// Schema, file hashes, job graph and portability (no absolute paths, globs or stamped names).
[[nodiscard]] json::Value validateDataPack(const std::filesystem::path& folder);

// Recomputes pack.sha256 and the manifest's file table (for pack authors).
void refreshDataPackHashes(const std::filesystem::path& folder);

[[nodiscard]] std::filesystem::path dataPackRootJob(const std::filesystem::path& folder);

// Compares a finished run's ledger with the pack's expected digests.
// verdict: verified_exact | numerically_equivalent | failed | non_reproducing | unsealed.
[[nodiscard]] json::Value verifyDataPackRun(const std::filesystem::path& folder,
                                            const std::filesystem::path& ledger);

// Copies the pack to `output`, recording the run's digests and this binary's build identity.
[[nodiscard]] json::Value sealDataPack(const std::filesystem::path& folder,
                                       const std::filesystem::path& ledger,
                                       const std::filesystem::path& output);

// Writes the compact zip DynSysLab-datapack-<name>-<date>-v<version>.zip into `directory`.
[[nodiscard]] json::Value exportDataPack(const std::filesystem::path& folder,
                                         const std::filesystem::path& directory);

// Walks an archive's recorded lineage into a draft pack description, listing every gap.
[[nodiscard]] json::Value traceArchiveLineage(const std::filesystem::path& archive);

}  // namespace dynsyslab::workflow
