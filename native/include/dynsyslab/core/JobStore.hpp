#pragma once
// Durable native job store (SQLite): submission, scheduling and job state.
// Scheduling follows the batch-scheduler conventions of Slurm (Yoo, Jette & Grondona 2003,
// "SLURM: Simple Linux Utility for Resource Management", JSSPP, LNCS 2862): dependency
// conditions afterok/afterany/afternotok/aftercorr, job arrays with a concurrency cap, and
// admission against declared thread and memory requests. See docs/architecture/JOB_QUEUE.md.

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace dynsyslab::core {

// queued -> running -> completed | failed | cancelled. paused jobs wait for `resume`;
// quarantined jobs (a pause the worker never acknowledged) wait for review; reserved
// jobs sit in the side tray until `release`.
enum class JobStatus { queued, running, completed, failed, cancelled, paused, quarantined, reserved };

enum class DependencyCondition { afterok, afterany, afternotok, aftercorr };

// For aftercorr, `job` names a job array: array task i waits for task i of that array.
struct JobDependency {
    std::string job;
    DependencyCondition condition{DependencyCondition::afterok};
};

// Declared requests used for admission. Zero means "not declared" (one thread, no memory).
struct JobResources {
    std::size_t threads{};
    std::size_t memoryBytes{};
};

struct JobSubmission {
    std::string label;
    std::string kind;
    std::string manifestJson;
    std::size_t retryLimit{};  // failed attempts requeued while attempts <= retryLimit
    std::vector<JobDependency> dependencies;
    JobResources resources;
    std::string arrayId;       // empty: not an array task
    std::size_t arrayIndex{};
    std::size_t arrayMaxConcurrent{};  // 0: no cap
    std::string specKey;       // spec claim: jobs with one key never run at the same time
    bool reserved{};           // side tray
};

struct JobRecord {
    std::string id;
    std::string label;
    std::string kind;
    std::string manifestJson;
    JobStatus status{JobStatus::queued};
    double progress{};
    std::string message;
    std::string worker;
    std::string createdUtc;
    std::string updatedUtc;
    std::size_t attempts{};
    std::size_t retryLimit{};
    std::vector<JobDependency> dependencies;
    JobResources resources;
    std::string arrayId;
    std::size_t arrayIndex{};
    std::size_t arrayMaxConcurrent{};
    std::string specKey;
    std::uint64_t processId{};  // worker process, for the pause watchdog
    std::string host;
    bool pauseRequested{};
};

// Pool limits checked before a job is claimed. Zero means unlimited.
struct AdmissionPolicy {
    std::size_t maxRunning{};
    std::size_t maxThreads{};
    std::size_t maxMemoryBytes{};
};

struct WatchdogFinding {
    std::string job;
    std::string classification;  // finishing_current_work | rogue_suspected | exited_without_acknowledgement
    bool workerAlive{};
    double heartbeatAgeSeconds{};
    bool quarantined{};
};

class JobStore {
public:
    explicit JobStore(const std::filesystem::path& path);
    ~JobStore();
    JobStore(const JobStore&) = delete;
    JobStore& operator=(const JobStore&) = delete;

    [[nodiscard]] std::string submit(const std::string& label, const std::string& kind,
                                     const std::string& manifestJson);
    [[nodiscard]] std::string submit(const JobSubmission& submission);
    // Claims the oldest queued job whose dependencies are met, whose array is under its cap,
    // whose spec key is not running elsewhere and which fits the pool. Jobs whose dependencies
    // can no longer be met are cancelled with the reason.
    [[nodiscard]] std::optional<JobRecord> claimNext(
        const std::string& worker, const AdmissionPolicy& policy = {},
        std::uint64_t processId = 0);
    void updateProgress(const std::string& id, double progress, const std::string& message);
    void complete(const std::string& id, const std::string& message);
    // Returns true when the failure was requeued under the job's retry limit.
    bool fail(const std::string& id, const std::string& message);
    void cancel(const std::string& id, const std::string& message = {});
    // queued/reserved -> paused at once; running -> pause requested at the next safe boundary.
    [[nodiscard]] JobStatus pause(const std::string& id);
    void acknowledgePause(const std::string& id, const std::string& message);
    void resume(const std::string& id);                 // paused or quarantined -> queued
    void release(const std::string& id);                // reserved -> queued
    void remove(const std::string& id);                 // deletes a reserved job
    void quarantine(const std::string& id, const std::string& message);
    [[nodiscard]] std::size_t recoverStale(double ageSeconds);
    // Classifies running jobs with a pause request; quarantines exited or silent workers.
    [[nodiscard]] std::vector<WatchdogFinding> watchPausedJobs(double graceSeconds);
    [[nodiscard]] std::vector<JobRecord> list(
        std::size_t limit = 100, std::optional<JobStatus> status = std::nullopt) const;
    [[nodiscard]] std::optional<JobRecord> find(const std::string& id) const;
    // Newest job with this spec key and status (completed by default).
    [[nodiscard]] std::optional<JobRecord> findBySpec(
        const std::string& specKey, JobStatus status = JobStatus::completed) const;
    [[nodiscard]] std::vector<JobRecord> arrayTasks(const std::string& arrayId) const;

    struct Implementation;

private:
    std::unique_ptr<Implementation> implementation_;
};

[[nodiscard]] std::string jobStatusName(JobStatus status);
[[nodiscard]] JobStatus parseJobStatus(const std::string& value);
[[nodiscard]] std::string dependencyConditionName(DependencyCondition condition);
[[nodiscard]] DependencyCondition parseDependencyCondition(const std::string& value);

// Scheduling rules, exposed for tests and planners.
enum class DependencyState { waiting, satisfied, impossible };
[[nodiscard]] DependencyState dependencyState(DependencyCondition condition,
                                              std::optional<JobStatus> status);
// Empty when the candidate fits the pool given the running jobs; otherwise the limit name.
[[nodiscard]] std::string admissionDenial(const JobResources& candidate,
                                          const std::vector<JobRecord>& running,
                                          const AdmissionPolicy& policy);
[[nodiscard]] bool processAlive(std::uint64_t processId);
[[nodiscard]] std::uint64_t currentProcessId();
[[nodiscard]] std::string currentHost();

}  // namespace dynsyslab::core
