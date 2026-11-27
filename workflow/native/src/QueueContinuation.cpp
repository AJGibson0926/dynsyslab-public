#include "QueueContinuation.hpp"

#include <dynsyslab/core/LongPath.hpp>
#include <dynsyslab/core/Text.hpp>

#include <chrono>
#include <fstream>
#include <stdexcept>
#include <utility>

namespace dynsyslab::workflow::detail {
namespace {

std::filesystem::path sidecarPath(const std::filesystem::path& root, const std::string& id) {
    return root / L"job_manager" / L"continuations" / core::wide(id + ".json");
}

void writeAtomic(const std::filesystem::path& path, const std::string& text) {
    std::filesystem::create_directories(path.parent_path());
    auto temporary = path;
    temporary += L".partial";
    {
        std::ofstream out(core::extendedLengthPath(temporary), std::ios::binary | std::ios::trunc);
        if (!out) throw std::runtime_error("unable to write the queue continuation record");
        out << text;
    }
    std::filesystem::rename(core::extendedLengthPath(temporary), core::extendedLengthPath(path));
}

// The solve's outcome record beside its hierarchy manifest, or null.
json::Value solveOutcome(const RunResult& run) {
    if (run.kind != "solve" || run.manifest.empty()) return {};
    auto path = run.manifest;
    path += L".outcome.json";
    if (!std::filesystem::is_regular_file(core::extendedLengthPath(path))) return {};
    return json::load(core::extendedLengthPath(path).wstring());
}

json::Value pathsJson(const std::vector<std::filesystem::path>& paths) {
    json::Value::Array result;
    for (const auto& path : paths) result.emplace_back(path.generic_u8string());
    return json::Value(std::move(result));
}

}  // namespace

bool checkpointPausable(const json::Value& job) {
    return job.at("kind").string() == "solve" && !job.at("campaign").at("enabled").boolean(false)
        && job.at("history").at("pause_at_checkpoint").boolean(false);
}

std::filesystem::path queuePauseFile(const std::filesystem::path& root, const std::string& id) {
    return root / L"job_manager" / L"pause" / core::wide(id + ".pause");
}

json::Value queuedJob(const json::Value& submitted, const std::filesystem::path& root,
                      const std::string& id, double wallTime) {
    if (submitted.at("kind").string() != "solve" || submitted.at("campaign").at("enabled").boolean(false))
        return submitted;
    json::Value job = submitted;
    if (checkpointPausable(job) && std::as_const(job).at("history").at("pause_file").isNull()) {
        const auto pause = queuePauseFile(root, id);
        std::error_code ignored;
        std::filesystem::remove(core::extendedLengthPath(pause), ignored);
        job.at("history").at("pause_file") = pause.generic_u8string();
    }
    const double margin = job.at("queue").at("checkpoint_before_deadline_seconds").number(0.0);
    if (margin < 0.0) throw std::invalid_argument("queue.checkpoint_before_deadline_seconds must be nonnegative");
    if (margin > 0.0 && wallTime > 0.0) {
        if (margin >= wallTime)
            throw std::invalid_argument("queue.checkpoint_before_deadline_seconds must be below wall_time_seconds");
        const double now = std::chrono::duration<double>(
            std::chrono::system_clock::now().time_since_epoch()).count();
        job.at("history").at("deadline_unix_seconds") = now + wallTime - margin;
    }
    return job;
}

CheckpointedRun runCheckpointed(const json::Value& job, const std::filesystem::path& root,
                                const std::string& id, const std::filesystem::path& progressPath,
                                const core::ExecutionContext* execution) {
    const auto sidecar = sidecarPath(root, id);
    json::Value::Array pending;
    std::vector<std::filesystem::path> completed;
    double rounds = 0.0;
    const bool resuming = std::filesystem::is_regular_file(core::extendedLengthPath(sidecar));
    if (resuming) {
        const auto state = json::load(core::extendedLengthPath(sidecar).wstring());
        pending = state.at("pending_jobs").array();
        for (const auto& path : state.at("completed_archives").array())
            completed.push_back(std::filesystem::u8path(path.string()));
        rounds = state.at("rounds").number(0.0);
    } else {
        pending.push_back(job);
    }
    const auto save = [&](const json::Value::Array& remaining) {
        writeAtomic(sidecar, json::dump(json::Value(json::Value::Object{
            {"schema", "dynsyslab.queue-continuation/v1"}, {"job_id", id}, {"rounds", rounds},
            {"completed_archives", pathsJson(completed)},
            {"pending_jobs", json::Value(remaining)}})));
    };
    CheckpointedRun result;
    RunResult last;
    for (std::size_t index = 0; index < pending.size(); ++index) {
        json::Value next = pending[index];
        // Continuation jobs carry this attempt's history keys (pause file, deadline, intervals).
        if (next.at("kind").string() == "solve" && !job.at("history").isNull())
            next.at("history") = job.at("history");
        last = executeJob(next, root, progressPath, execution);
        completed.insert(completed.end(), last.archives.begin(), last.archives.end());
        json::Value::Array remaining;
        const auto outcome = solveOutcome(last);
        const bool paused = outcome.at("paused").boolean(false);
        if (paused)
            for (const auto& resume : outcome.at("resume_jobs").array()) remaining.push_back(resume);
        for (std::size_t later = index + 1; later < pending.size(); ++later) remaining.push_back(pending[later]);
        if (paused) {
            rounds += 1.0;
            save(remaining);
            result.paused = true;
            break;
        }
        if (!remaining.empty()) save(remaining);
    }
    if (!result.paused) {
        std::error_code ignored;
        std::filesystem::remove(core::extendedLengthPath(sidecar), ignored);
    }
    result.run = std::move(last);
    result.run.archives = completed;
    result.record = json::Value::Object{{"schema", "dynsyslab.queue-continuation/v1"},
        {"rounds", rounds}, {"resumed", resuming}, {"paused", result.paused}};
    return result;
}

}  // namespace dynsyslab::workflow::detail
