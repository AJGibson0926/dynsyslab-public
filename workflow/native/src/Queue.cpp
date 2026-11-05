#include <dynsyslab/workflow/Queue.hpp>

#include <dynsyslab/core/Text.hpp>
#include <dynsyslab/workflow/Plan.hpp>

#include <chrono>
#include <stdexcept>
#include <utility>

namespace dynsyslab::workflow {
namespace {

json::Value paths(const std::vector<std::filesystem::path>& values) {
    json::Value::Array result;
    for (const auto& value : values) result.emplace_back(value.generic_u8string());
    return json::Value(std::move(result));
}

bool wasCancelled(dynsyslab::core::JobStore& store, const std::string& id) {
    const auto current = store.find(id);
    return !current || current->status == dynsyslab::core::JobStatus::cancelled;
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

JobQueue::JobQueue(std::filesystem::path archiveRoot)
    : archiveRoot_(std::move(archiveRoot)),
      database_(archiveRoot_ / L"job_manager" / L"state.sqlite") {}

std::string JobQueue::submit(const json::Value& job) {
    validateJob(job);
    (void)planJob(job);
    dynsyslab::core::JobStore store(database_);
    return store.submit(job.at("label").string("job"),
                        job.at("kind").string(), json::dump(job));
}

std::vector<dynsyslab::core::JobRecord> JobQueue::list(std::size_t limit) const {
    dynsyslab::core::JobStore store(database_);
    return store.list(limit);
}

std::optional<dynsyslab::core::JobRecord> JobQueue::find(
    const std::string& id) const {
    dynsyslab::core::JobStore store(database_);
    return store.find(id);
}

dynsyslab::core::JobRecord JobQueue::cancel(const std::string& id) {
    dynsyslab::core::JobStore store(database_);
    const auto current = store.find(id);
    if (!current) throw std::invalid_argument("queued job does not exist: " + id);
    if (current->status != dynsyslab::core::JobStatus::queued
        && current->status != dynsyslab::core::JobStatus::running)
        throw std::invalid_argument("only queued or running jobs may be cancelled");
    store.cancel(id, "cancelled by native UI or CLI");
    return *store.find(id);
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
    dynsyslab::core::CancellationProbe externalCancellation) {
    dynsyslab::core::JobStore store(database_);
    const auto claimed = store.claimNext(worker);
    if (!claimed) return {};

    bool cancelled = false;
    auto lastPoll = std::chrono::steady_clock::time_point{};
    dynsyslab::core::ExecutionContext execution;
    execution.limits = limits;
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
        if (lastPoll.time_since_epoch().count() == 0
            || now - lastPoll >= std::chrono::milliseconds(100)) {
            cancelled = wasCancelled(store, claimed->id);
            lastPoll = now;
        }
        return cancelled;
    };
    execution.progress = [&](const dynsyslab::core::ExecutionProgress& event) {
        store.updateProgress(claimed->id, event.fraction, event.message);
    };

    try {
        store.updateProgress(claimed->id, 0.0, "native execution started");
        const auto activeProgress = progressPath.empty()
            ? archiveRoot_ / L"job_manager" / L"progress"
                / dynsyslab::core::wide(claimed->id + ".json")
            : progressPath;
        auto run = executeJob(json::parse(claimed->manifestJson), archiveRoot_,
                              activeProgress, &execution);
        if (wasCancelled(store, claimed->id))
            return {QueueRunResult::State::cancelled, store.find(claimed->id), {}};
        store.complete(claimed->id, json::dump(queueRunJson(run)));
        return {QueueRunResult::State::completed, store.find(claimed->id),
                std::move(run)};
    } catch (const dynsyslab::core::ExecutionCancelled&) {
        if (!wasCancelled(store, claimed->id))
            store.cancel(claimed->id, "native execution cancelled");
        return {QueueRunResult::State::cancelled, store.find(claimed->id), {}};
    } catch (const std::exception& error) {
        if (!wasCancelled(store, claimed->id)) store.fail(claimed->id, error.what());
        throw;
    }
}

}  // namespace dynsyslab::workflow
