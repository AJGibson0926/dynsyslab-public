#pragma once

#include <dynsyslab/workflow/Json.hpp>

namespace dynsyslab::workflow {

// Validate a job against the compiled registries and return resource/layout
// estimates without running numerical work or writing archives.
[[nodiscard]] json::Value planJob(const json::Value& job);

}  // namespace dynsyslab::workflow
