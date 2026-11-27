#pragma once
// Verified input staging for a queue job (port of job_manager/staging.py): the files named in
// `queue.staging.inputs` are copied into job_manager/staging/<job id>/, each copy's SHA-256 is
// checked against its source, and staging-manifest.json (dynsyslab.staging-manifest/v1) lists
// source, staged path, digest and size. Primary archives are only read. Remote commands get
// the folder as {staging_dir}; `cleanup_scratch` (default true) removes it afterwards.

#include <dynsyslab/workflow/Json.hpp>

#include <filesystem>
#include <string>

namespace dynsyslab::workflow::detail {

// Null when the job declares no staging inputs.
[[nodiscard]] json::Value stageInputs(const json::Value& job, const std::filesystem::path& archiveRoot,
                                      const std::string& jobId);
void cleanupStaging(const json::Value& manifest);

}  // namespace dynsyslab::workflow::detail
