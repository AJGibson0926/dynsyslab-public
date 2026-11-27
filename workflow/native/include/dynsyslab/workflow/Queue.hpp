#pragma once
// Native job queue: submission (retry, dependencies, arrays, spec claims, side tray, reuse),
// execution of the next admissible job, the resident `serve` loop and command-template
// remote execution. Job keys and CLI flags: docs/architecture/JOB_QUEUE.md.

#include <dynsyslab/core/Execution.hpp>
#include <dynsyslab/core/JobStore.hpp>
#include <dynsyslab/workflow/Json.hpp>
#include <dynsyslab/workflow/Run.hpp>

#include <cstddef>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace dynsyslab::workflow {

struct QueueRunResult {
    // requeued: the job failed and went back to the queue under its retry limit.
    // paused: the job stopped at a safe boundary after `queue pause`.
    enum class State { empty, completed, cancelled, requeued, paused };
    State state{State::empty};
    std::optional<dynsyslab::core::JobRecord> job;
    std::optional<RunResult> run;
    std::string error;
};

struct QueueSubmitResult {
    std::vector<std::string> jobIds;
    std::string arrayId;
    std::string campaignId;
    std::string reusedJob;            // an identical completed job; nothing new was queued
    json::Value extensionCandidates;  // shorter or longer completed runs of the same solve
    std::string extendedFrom;         // queue.reuse "extend": the shorter job continued
    [[nodiscard]] json::Value toJson() const;
};

struct QueueServeOptions {
    std::string worker{"dynsyslab-queue-serve"};
    std::size_t workers{1};
    dynsyslab::core::AdmissionPolicy admission;
    dynsyslab::core::ExecutionLimits limits;
    double pollSeconds{2.0};
    std::size_t maxJobs{};         // 0: unlimited
    double idleExitSeconds{-1.0};  // < 0: never exit when idle
    double watchdogGraceSeconds{120.0};
    // When set, each claimed job is written to job_manager/remote/<id>.json and the command
    // runs with {job_file}, {job_id}, {archive_root} and {result_file} substituted.
    std::string commandTemplate;
};

struct QueueServeSummary {
    std::size_t completed{}, failed{}, requeued{}, cancelled{}, paused{}, quarantined{};
    [[nodiscard]] json::Value toJson() const;
};

class JobQueue {
public:
    explicit JobQueue(std::filesystem::path archiveRoot);

    // Returns the first queued id (or the reused id); see submitDetailed for arrays.
    [[nodiscard]] std::string submit(const json::Value& job);
    [[nodiscard]] QueueSubmitResult submitDetailed(const json::Value& job);
    [[nodiscard]] std::vector<dynsyslab::core::JobRecord> list(
        std::size_t limit = 100) const;
    [[nodiscard]] std::optional<dynsyslab::core::JobRecord> find(
        const std::string& id) const;
    [[nodiscard]] dynsyslab::core::JobRecord cancel(const std::string& id);
    [[nodiscard]] dynsyslab::core::JobRecord pause(const std::string& id);
    [[nodiscard]] dynsyslab::core::JobRecord resume(const std::string& id);
    [[nodiscard]] dynsyslab::core::JobRecord release(const std::string& id);
    void remove(const std::string& id);
    [[nodiscard]] std::vector<dynsyslab::core::JobRecord> tray() const;
    [[nodiscard]] std::vector<dynsyslab::core::WatchdogFinding> watchdog(double graceSeconds);
    [[nodiscard]] std::size_t recoverStale(double ageSeconds);
    // `limits` caps the job's own execution limits (job key `execution`); `admission` bounds
    // what may run next to the jobs already running.
    [[nodiscard]] QueueRunResult runNext(
        const std::string& worker,
        const std::filesystem::path& progressPath = {},
        dynsyslab::core::ExecutionLimits limits = {},
        dynsyslab::core::CancellationProbe externalCancellation = {},
        const dynsyslab::core::AdmissionPolicy& admission = {});
    [[nodiscard]] QueueServeSummary serve(const QueueServeOptions& options,
                                          dynsyslab::core::CancellationProbe stop = {});

    [[nodiscard]] const std::filesystem::path& archiveRoot() const { return archiveRoot_; }
    [[nodiscard]] const std::filesystem::path& database() const { return database_; }

private:
    QueueRunResult runRemote(const dynsyslab::core::JobRecord& job,
                             const std::string& commandTemplate);
    std::filesystem::path archiveRoot_;
    std::filesystem::path database_;
};

[[nodiscard]] json::Value queueRunJson(const RunResult& run);
[[nodiscard]] json::Value queueJobJson(const dynsyslab::core::JobRecord& job);
// Smaller of the two limits per field, treating zero as unlimited.
[[nodiscard]] dynsyslab::core::ExecutionLimits combinedLimits(
    dynsyslab::core::ExecutionLimits first, dynsyslab::core::ExecutionLimits second);
// Limits declared by a job: `execution.max_threads`, `execution.max_memory_mib`.
[[nodiscard]] dynsyslab::core::ExecutionLimits jobExecutionLimits(const json::Value& job);
// Stable key of a job's scientific content (label, queue options and ids removed).
[[nodiscard]] std::string jobSpecKey(const json::Value& job);

// Campaign registry reports (registry rows are written by campaigns and array submission).
[[nodiscard]] json::Value campaignList(const std::filesystem::path& archiveRoot);
[[nodiscard]] json::Value campaignStatus(const std::filesystem::path& archiveRoot,
                                         const std::string& campaignId);
// Checks every point of the original plan against the published chunk manifests: each
// point covered once, with the planned coordinates and existing solution archives.
[[nodiscard]] json::Value campaignValidation(const std::filesystem::path& archiveRoot,
                                             const std::string& campaignId);
// Resubmits, from the original plan, every chunk that has no published manifest and whose
// queue job failed, was cancelled, was quarantined or is missing. dynsyslab.campaign-repair/v1.
[[nodiscard]] json::Value campaignRepair(const std::filesystem::path& archiveRoot,
                                         const std::string& campaignId);
// Splits every unstarted (queued, reserved or paused) unpublished chunk larger than
// pointsPerChunk into smaller queue jobs over the same points, cancels the original job and
// records the split in job_manager/restructure/ (dynsyslab.job-restructure/v1). The plan and
// its point ordering are unchanged. Port of JobManager.restructure.
[[nodiscard]] json::Value campaignRestructure(const std::filesystem::path& archiveRoot,
    const std::string& campaignId, std::size_t pointsPerChunk, const std::string& reason = {});
// Read-only coverage review before (or while) a campaign runs: every planned point covered by
// exactly one chunk, every chunk backed by a queue job or a manifest that exists, and the
// chunks needing `campaign repair`. Written to job_manager/reviews/
// (dynsyslab.campaign-review/v1). Port of review_parallel_transition.
[[nodiscard]] json::Value campaignReview(const std::filesystem::path& archiveRoot,
                                         const std::string& campaignId);
// Local resources for parallel work without running a solve: logical cores, free archive
// space, running jobs and a recommended worker count. dynsyslab.queue-resources/v1.
[[nodiscard]] json::Value queueResources(const std::filesystem::path& archiveRoot);

}  // namespace dynsyslab::workflow
