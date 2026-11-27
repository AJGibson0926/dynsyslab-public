#pragma once
// Campaign design for solve jobs: quality, execution and checkpoint profiles; linspace and
// driver.* axes; points_range chunks; and cycle-based time specifications.

#include <dynsyslab/workflow/Json.hpp>

#include <cstddef>
#include <string>

namespace dynsyslab::workflow::detail {

// Returns the job with its campaign design expanded: every axis carries explicit `values`
// and a native path, profile defaults are filled in, and a cycle-based time without a swept
// frequency is resolved to start/stop/output_step. Jobs without these keys are unchanged.
[[nodiscard]] json::Value resolveJobDesign(const json::Value& job);

// Resolves `time.transient_cycles`, `time.retained_cycles` (or `time.cycle_profile`) and
// `time.samples_per_cycle` against the driver frequency into start/stop/output_step, and
// records `time.retained_start`. Idempotent; a job without cycle keys is returned unchanged.
[[nodiscard]] json::Value resolveCycleTime(const json::Value& job);

// Campaign points handled by one queue task: `systems_per_job` from the job or its profile.
[[nodiscard]] std::size_t campaignSystemsPerJob(const json::Value& job);

// The named profile tables (dynsyslab.campaign-profiles/v1); public as workflow::campaignProfiles.
[[nodiscard]] json::Value campaignProfiles();

}  // namespace dynsyslab::workflow::detail
