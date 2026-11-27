#pragma once
// Declarative history (checkpoint) policy, decided before execution from the job's runtime
// estimate (the native port of job_manager/history_policy.py). The solver never decides
// whether a run is worth checkpointing.
//
// Policy, from `history_policy` (any solve job) or `campaign.checkpoint_policy` (set by
// `campaign.checkpoint_profile`):
//   enabled, estimated_seconds (default queue.estimated_seconds, then queue.wall_time_seconds),
//   if_estimated_seconds_over, when_estimate_unknown ("disable" default | "enable"),
//   mode ("checkpoint" default | "telemetry": progress only, no checkpoints),
//   checkpoint_interval (model time), checkpoint_interval_seconds or interval_seconds
//   (wall clock), pause_at_checkpoint.
// The resolved settings go to job.history (checkpoint_interval, checkpoint_wall_seconds,
// pause_at_checkpoint) and the decision to job.history.policy_resolution
// {enabled, reason, estimated_seconds, threshold_seconds, mode}. Idempotent.

#include <dynsyslab/workflow/Json.hpp>

namespace dynsyslab::workflow::detail {

[[nodiscard]] json::Value applyHistoryPolicy(const json::Value& job);

}  // namespace dynsyslab::workflow::detail
