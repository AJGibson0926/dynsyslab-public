#pragma once
// Queue-side solver checkpoints: the deadline policy (checkpoint and pause before
// queue.wall_time_seconds, requeue, continue from the checkpoint) and pause at a checkpoint
// (history.pause_at_checkpoint). A paused solve's continuation jobs (the outcome record's
// resume_jobs, SolveContinuation.hpp) are kept in job_manager/continuations/<id>.json; the
// next attempt of the same queue job runs them instead of restarting, so the job id is stable.
// Port of job_manager/scheduling.py DeadlinePolicy and the pause-at-checkpoint path.

#include <dynsyslab/core/Execution.hpp>
#include <dynsyslab/workflow/Json.hpp>
#include <dynsyslab/workflow/Run.hpp>

#include <filesystem>
#include <string>

namespace dynsyslab::workflow::detail {

// A non-campaign solve with history.pause_at_checkpoint: a queue pause goes through the
// solve's pause file so the solve stops at its next checkpoint instead of being cancelled.
[[nodiscard]] bool checkpointPausable(const json::Value& job);

[[nodiscard]] std::filesystem::path queuePauseFile(const std::filesystem::path& archiveRoot,
                                                   const std::string& jobId);

// The job as the queue runs it: history.pause_file for checkpoint-pausable jobs (a stale
// file from an earlier pause is removed) and, with queue.checkpoint_before_deadline_seconds
// and a wall time, history.deadline_unix_seconds = now + wall time - margin.
[[nodiscard]] json::Value queuedJob(const json::Value& job, const std::filesystem::path& archiveRoot,
                                    const std::string& jobId, double wallTimeSeconds);

struct CheckpointedRun {
    RunResult run;          // every archive of this job so far, earlier attempts included
    bool paused{};          // continuation jobs remain in the sidecar
    json::Value record;     // dynsyslab.queue-continuation/v1 summary for the job result
};

// Runs the job, or the continuation jobs a paused earlier attempt left.
[[nodiscard]] CheckpointedRun runCheckpointed(const json::Value& job,
    const std::filesystem::path& archiveRoot, const std::string& jobId,
    const std::filesystem::path& progressPath, const core::ExecutionContext* execution);

}  // namespace dynsyslab::workflow::detail
