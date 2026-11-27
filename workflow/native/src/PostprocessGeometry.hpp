#pragma once
// Post-processing plan geometry and the dense memory preflight, computed from the source
// archives without running any method: samples, dimension and uniform clock of each source;
// Hankel shape and retained-rank bound; SINDy feature count; a dense working-set estimate.

#include <dynsyslab/core/Execution.hpp>
#include <dynsyslab/workflow/Json.hpp>

namespace dynsyslab::workflow::detail {

// dynsyslab.postprocess-geometry/v1 with `valid`, `issues` and `estimated_peak_bytes`.
// Sources that do not exist yet (experiment placeholders) are listed as unresolved.
[[nodiscard]] json::Value postprocessGeometry(const json::Value& job);

// Refuses the job before any allocation when the estimate exceeds the memory budget: the
// smaller of `limits.maximumMemoryBytes`, `execution.max_memory_mib` and each stage's
// `settings.max_memory_mb`. Throws core::ExecutionLimitExceeded.
void requirePostprocessMemory(const json::Value& job, const dynsyslab::core::ExecutionLimits& limits);

}  // namespace dynsyslab::workflow::detail
