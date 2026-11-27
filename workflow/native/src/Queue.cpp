#include <dynsyslab/workflow/Queue.hpp>

#include <dynsyslab/core/Text.hpp>
#include <dynsyslab/workflow/Plan.hpp>

#include "QueueContinuation.hpp"
#include "QueueStaging.hpp"

#include <chrono>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <utility>

namespace dynsyslab::workflow {
namespace {

using dynsyslab::core::JobStatus;

json::Value paths(const std::vector<std::filesystem::path>& values) {
    json::Value::Array result;
    for (const auto& value : values) result.emplace_back(value.generic_u8string());
    return json::Value(std::move(result));
}

bool wasCancelled(dynsyslab::core::JobStore& store, const std::string& id) {
    const auto current = store.find(id);
    return !current || current->status == JobStatus::cancelled
        || current->status == JobStatus::quarantined;
}

std::size_t smaller(std::size_t first, std::size_t second) {
    if (!first) return second;
    if (!second) return first;
    return std::min(first, second);
}

}  // namespace

json::Value queueRunJson(const RunResult& run) {
    json::Value::Object result;
    result["schema"] = "dynsyslab.workflow-result/v1";
    result["kind"] = run.kind;
    result["archives"] = paths(run.archives);
    result["summary"] = dynsyslab::core::narrow(run.summary);
    if (!run.manifest.empty())
        result[run.kind == "campaign" || run.kind == "krom_control_search"
            || run.kind == "solve_set"
            ? "manifest" : "hierarchy_manifest"] = run.manifest.generic_u8string();
    return json::Value(std::move(result));
}

json::Value queueJobJson(const dynsyslab::core::JobRecord& job) {
    json::Value::Object value;
    value["id"] = job.id; value["label"] = job.label; value["kind"] = job.kind;
    value["status"] = dynsyslab::core::jobStatusName(job.status);
    value["progress"] = job.progress; value["message"] = job.message;
    value["worker"] = job.worker; value["created_utc"] = job.createdUtc;
    value["updated_utc"] = job.updatedUtc;
    value["attempts"] = static_cast<double>(job.attempts);
    value["retry_limit"] = static_cast<double>(job.retryLimit);
    json::Value::Array dependencies;
    for (const auto& item : job.dependencies)
        dependencies.emplace_back(json::Value::Object{{"job", item.job},
            {"condition", dynsyslab::core::dependencyConditionName(item.condition)}});
    value["dependencies"] = json::Value(std::move(dependencies));
    value["resources"] = json::Value::Object{
        {"threads", static_cast<double>(job.resources.threads)},
        {"memory_mib", static_cast<double>(job.resources.memoryBytes) / (1024.0 * 1024.0)}};
    if (!job.arrayId.empty())
        value["array"] = json::Value::Object{{"id", job.arrayId},
            {"index", static_cast<double>(job.arrayIndex)},
            {"max_concurrent", static_cast<double>(job.arrayMaxConcurrent)}};
    if (!job.specKey.empty()) value["spec_key"] = job.specKey;
    if (job.pauseRequested) value["pause_requested"] = true;
    return json::Value(std::move(value));
}

dynsyslab::core::ExecutionLimits combinedLimits(dynsyslab::core::ExecutionLimits first,
                                                dynsyslab::core::ExecutionLimits second) {
    return {smaller(first.maximumThreads, second.maximumThreads),
            smaller(first.maximumMemoryBytes, second.maximumMemoryBytes)};
}

dynsyslab::core::ExecutionLimits jobExecutionLimits(const json::Value& job) {
    const auto& execution = job.at("execution");
    const double threads = execution.at("max_threads").number(0.0);
    const double memory = execution.at("max_memory_mib").number(0.0);
    if (threads < 0.0 || threads != std::floor(threads) || memory < 0.0)
        throw std::invalid_argument(
            "execution.max_threads must be a whole number and execution.max_memory_mib nonnegative");
    return {static_cast<std::size_t>(threads),
            static_cast<std::size_t>(memory * 1024.0 * 1024.0)};
}

JobQueue::JobQueue(std::filesystem::path archiveRoot)
    : archiveRoot_(std::move(archiveRoot)),
      database_(archiveRoot_ / L"job_manager" / L"state.sqlite") {}

std::string JobQueue::submit(const json::Value& job) {
    const auto result = submitDetailed(job);
    return result.reusedJob.empty() ? result.jobIds.front() : result.reusedJob;
}

std::vector<dynsyslab::core::JobRecord> JobQueue::list(std::size_t limit) const {
    dynsyslab::core::JobStore store(database_);
    return store.list(limit);
}

std::optional<dynsyslab::core::JobRecord> JobQueue::find(const std::string& id) const {
    dynsyslab::core::JobStore store(database_);
    return store.find(id);
}

dynsyslab::core::JobRecord JobQueue::cancel(const std::string& id) {
    dynsyslab::core::JobStore store(database_);
    const auto current = store.find(id);
    if (!current) throw std::invalid_argument("queued job does not exist: " + id);
    if (current->status != JobStatus::queued && current->status != JobStatus::running
        && current->status != JobStatus::paused && current->status != JobStatus::reserved
        && current->status != JobStatus::quarantined)
        throw std::invalid_argument("only unfinished jobs may be cancelled");
    store.cancel(id, "cancelled by native UI or CLI");
    return *store.find(id);
}

dynsyslab::core::JobRecord JobQueue::pause(const std::string& id) {
    dynsyslab::core::JobStore store(database_);
    (void)store.pause(id);
    return *store.find(id);
}

dynsyslab::core::JobRecord JobQueue::resume(const std::string& id) {
    dynsyslab::core::JobStore store(database_);
    store.resume(id);
    return *store.find(id);
}

dynsyslab::core::JobRecord JobQueue::release(const std::string& id) {
    dynsyslab::core::JobStore store(database_);
    store.release(id);
    return *store.find(id);
}

void JobQueue::remove(const std::string& id) {
    dynsyslab::core::JobStore store(database_);
    store.remove(id);
}

std::vector<dynsyslab::core::JobRecord> JobQueue::tray() const {
    dynsyslab::core::JobStore store(database_);
    return store.list(1000, JobStatus::reserved);
}

std::vector<dynsyslab::core::WatchdogFinding> JobQueue::watchdog(double graceSeconds) {
    if (graceSeconds < 0.0) throw std::invalid_argument("watchdog grace must be nonnegative");
    dynsyslab::core::JobStore store(database_);
    return store.watchPausedJobs(graceSeconds);
}

std::size_t JobQueue::recoverStale(double ageSeconds) {
    if (ageSeconds < 0.0)
        throw std::invalid_argument("stale recovery age must be nonnegative");
    dynsyslab::core::JobStore store(database_);
    return store.recoverStale(ageSeconds);
}

QueueRunResult JobQueue::runNext(const std::string& worker,
    const std::filesystem::path& progressPath,
    dynsyslab::core::ExecutionLimits limits,
    dynsyslab::core::CancellationProbe externalCancellation,
    const dynsyslab::core::AdmissionPolicy& admission) {
    dynsyslab::core::JobStore store(database_);
    const auto claimed = store.claimNext(worker, admission, dynsyslab::core::currentProcessId());
    if (!claimed) return {};

    bool cancelled = false, pauseSeen = false, timedOut = false;
    auto lastPoll = std::chrono::steady_clock::time_point{};
    // queue.wall_time_seconds: past it the job stops at its next checkpoint and fails, so it
    // is requeued under its retry limit (requeue on timeout; retry_limit = maximum requeues).
    // queue.checkpoint_before_deadline_seconds: a solve checkpoints and pauses that long before
    // it, and is requeued to continue from the checkpoint (QueueContinuation.hpp).
    double wallTime = 0.0;
    bool pauseThroughCheckpoint = false;
    try {
        const auto manifest = json::parse(claimed->manifestJson);
        wallTime = manifest.at("queue").at("wall_time_seconds").number(0.0);
        pauseThroughCheckpoint = detail::checkpointPausable(manifest);
    } catch (const std::exception&) { /* an unreadable manifest fails below, inside the try */ }
    const auto started = std::chrono::steady_clock::now();
    dynsyslab::core::ExecutionContext execution;
    execution.cancellation = [&] {
        if (externalCancellation && externalCancellation()) {
            if (!cancelled) {
                try { store.cancel(claimed->id, "cancelled from native workbench"); }
                catch (...) { /* The job may already have reached a terminal state. */ }
            }
            cancelled = true;
            return true;
        }
        const auto now = std::chrono::steady_clock::now();
        if (wallTime > 0.0 && std::chrono::duration<double>(now - started).count() > wallTime) {
            timedOut = true;
            return true;
        }
        if (lastPoll.time_since_epoch().count() == 0
            || now - lastPoll >= std::chrono::milliseconds(100)) {
            const auto current = store.find(claimed->id);
            cancelled = !current || current->status == JobStatus::cancelled
                || current->status == JobStatus::quarantined;
            pauseSeen = current && current->pauseRequested;
            lastPoll = now;
            if (pauseSeen && pauseThroughCheckpoint && !cancelled) {
                // The solve stops at its next checkpoint and keeps its samples.
                const auto file = detail::queuePauseFile(archiveRoot_, claimed->id);
                std::filesystem::create_directories(file.parent_path());
                std::ofstream touch(file);
                return false;
            }
        }
        if (pauseSeen && pauseThroughCheckpoint && !cancelled) return false;
        return cancelled || pauseSeen;
    };
    execution.progress = [&](const dynsyslab::core::ExecutionProgress& event) {
        store.updateProgress(claimed->id, event.fraction, event.message);
    };

    try {
        const auto job = detail::queuedJob(json::parse(claimed->manifestJson), archiveRoot_,
                                           claimed->id, wallTime);
        execution.limits = combinedLimits(limits, jobExecutionLimits(job));
        const auto previousThreads = dynsyslab::core::processThreadLimit();
        if (execution.limits.maximumThreads)
            dynsyslab::core::setProcessThreadLimit(execution.limits.maximumThreads);
        struct RestoreThreads {
            std::size_t value;
            ~RestoreThreads() { dynsyslab::core::setProcessThreadLimit(value); }
        } restore{previousThreads};
        store.updateProgress(claimed->id, 0.0, "native execution started");
        const auto activeProgress = progressPath.empty()
            ? archiveRoot_ / L"job_manager" / L"progress"
                / dynsyslab::core::wide(claimed->id + ".json")
            : progressPath;
        const auto staging = detail::stageInputs(job, archiveRoot_, claimed->id);
        struct Cleanup { const json::Value& manifest; ~Cleanup() { detail::cleanupStaging(manifest); } } cleanup{staging};
        auto checkpointed = detail::runCheckpointed(job, archiveRoot_, claimed->id, activeProgress, &execution);
        auto& run = checkpointed.run;
        if (wasCancelled(store, claimed->id))
            return {QueueRunResult::State::cancelled, store.find(claimed->id), {}, {}};
        if (checkpointed.paused) {
            const double deadline = job.at("history").at("deadline_unix_seconds").number(0.0);
            const double now = std::chrono::duration<double>(
                std::chrono::system_clock::now().time_since_epoch()).count();
            if (deadline > 0.0 && now >= deadline && !pauseSeen) {
                const std::string message = "checkpointed before the wall-time deadline; continues from the checkpoint";
                if (store.fail(claimed->id, message))
                    return {QueueRunResult::State::requeued, store.find(claimed->id), {}, message};
                return {QueueRunResult::State::completed, store.find(claimed->id), {},
                        message + " (requeue limit reached)"};
            }
            store.acknowledgePause(claimed->id, "paused at a checkpoint; resumes from it");
            return {QueueRunResult::State::paused, store.find(claimed->id), {}, {}};
        }
        auto summary = queueRunJson(run);
        if (checkpointed.record.at("rounds").number(0.0) > 0.0) summary.at("continuation") = checkpointed.record;
        if (!staging.isNull()) summary.at("staging") = staging;
        store.complete(claimed->id, json::dump(summary));
        return {QueueRunResult::State::completed, store.find(claimed->id), std::move(run), {}};
    } catch (const dynsyslab::core::ExecutionCancelled&) {
        if (timedOut && !cancelled) {
            const std::string error = "wall time of " + std::to_string(wallTime) + " s exceeded";
            if (store.fail(claimed->id, error))
                return {QueueRunResult::State::requeued, store.find(claimed->id), {}, error};
            throw std::runtime_error(error);
        }
        if (pauseSeen && !cancelled) {
            // Not checkpoint-pausable: restarts from the beginning (or its continuation) on resume.
            store.acknowledgePause(claimed->id, "paused at a safe boundary");
            return {QueueRunResult::State::paused, store.find(claimed->id), {}, {}};
        }
        if (!wasCancelled(store, claimed->id))
            store.cancel(claimed->id, "native execution cancelled");
        return {QueueRunResult::State::cancelled, store.find(claimed->id), {}, {}};
    } catch (const std::exception& error) {
        if (wasCancelled(store, claimed->id)) throw;
        if (store.fail(claimed->id, error.what()))
            return {QueueRunResult::State::requeued, store.find(claimed->id), {}, error.what()};
        throw;
    }
}

}  // namespace dynsyslab::workflow
