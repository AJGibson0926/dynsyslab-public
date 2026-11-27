// Resident queue loop (`queue serve`) and command-template execution of claimed jobs, the
// native replacement for the Python daemon and its SSH/scheduler backends.
#include <dynsyslab/workflow/Queue.hpp>

#include <dynsyslab/core/Text.hpp>

#include "QueueStaging.hpp"

#include <atomic>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <thread>

namespace dynsyslab::workflow {
namespace {

using dynsyslab::core::JobStatus;

std::string replaceAll(std::string text, const std::string& token, const std::string& value) {
    for (std::size_t at = text.find(token); at != std::string::npos;
         at = text.find(token, at + value.size()))
        text.replace(at, token.size(), value);
    return text;
}

std::string readFile(const std::filesystem::path& path) {
    std::ifstream stream(path, std::ios::binary);
    std::ostringstream out; out << stream.rdbuf();
    return out.str();
}

}  // namespace

json::Value QueueServeSummary::toJson() const {
    return json::Value::Object{{"completed", static_cast<double>(completed)},
        {"failed", static_cast<double>(failed)}, {"requeued", static_cast<double>(requeued)},
        {"cancelled", static_cast<double>(cancelled)}, {"paused", static_cast<double>(paused)},
        {"quarantined", static_cast<double>(quarantined)}};
}

QueueRunResult JobQueue::runRemote(const dynsyslab::core::JobRecord& job,
                                   const std::string& commandTemplate) {
    dynsyslab::core::JobStore store(database_);
    const auto directory = archiveRoot_ / L"job_manager" / L"remote";
    std::filesystem::create_directories(directory);
    const auto jobFile = directory / dynsyslab::core::wide(job.id + ".json");
    const auto resultFile = directory / dynsyslab::core::wide(job.id + ".result.json");
    { std::ofstream out(jobFile, std::ios::binary | std::ios::trunc); out << job.manifestJson; }
    std::error_code ignored;
    std::filesystem::remove(resultFile, ignored);
    // queue.staging: verified copies of the inputs, removed after the command (QueueStaging.hpp).
    json::Value staging;
    try { staging = detail::stageInputs(json::parse(job.manifestJson), archiveRoot_, job.id); }
    catch (const std::exception& error) {
        const std::string message = std::string("input staging failed: ") + error.what();
        const bool requeued = store.fail(job.id, message);
        return {requeued ? QueueRunResult::State::requeued : QueueRunResult::State::completed,
                store.find(job.id), {}, message};
    }
    struct Cleanup { const json::Value& manifest; ~Cleanup() { detail::cleanupStaging(manifest); } } cleanup{staging};
    auto command = replaceAll(commandTemplate, "{job_file}", jobFile.u8string());
    command = replaceAll(command, "{staging_dir}",
                         staging.isNull() ? std::string{} : staging.at("staging_directory").string());
    command = replaceAll(command, "{job_id}", job.id);
    command = replaceAll(command, "{archive_root}", archiveRoot_.u8string());
    command = replaceAll(command, "{result_file}", resultFile.u8string());
    store.updateProgress(job.id, 0.0, "submitted through command template");
    const int status = std::system(command.c_str());
    const auto current = store.find(job.id);
    if (current && current->status == JobStatus::cancelled)
        return {QueueRunResult::State::cancelled, current, {}, {}};
    if (status == 0) {
        // The command may leave the remote workflow result in {result_file}.
        auto text = std::filesystem::is_regular_file(resultFile) ? readFile(resultFile)
            : std::string("{\"schema\":\"dynsyslab.workflow-result/v1\",\"remote\":true}");
        if (!staging.isNull()) {
            try { auto parsed = json::parse(text); parsed.at("staging") = staging; text = json::dump(parsed); }
            catch (const std::exception&) { /* a non-JSON remote result is kept as written */ }
        }
        store.complete(job.id, text);
        return {QueueRunResult::State::completed, store.find(job.id), {}, {}};
    }
    const std::string error = "command template exited with status " + std::to_string(status);
    const bool requeued = store.fail(job.id, error);
    return {requeued ? QueueRunResult::State::requeued : QueueRunResult::State::completed,
            store.find(job.id), {}, error};
}

QueueServeSummary JobQueue::serve(const QueueServeOptions& options,
                                  dynsyslab::core::CancellationProbe stop) {
    if (!options.workers) throw std::invalid_argument("queue serve needs at least one worker");
    if (options.pollSeconds < 0.0) throw std::invalid_argument("poll seconds must be nonnegative");
    QueueServeSummary summary;
    std::mutex lock;
    std::atomic<std::size_t> started{0};
    std::atomic<bool> finished{false};
    auto lastWork = std::chrono::steady_clock::now();
    const auto stopping = [&] {
        return finished.load() || (stop && stop())
            || (options.maxJobs && started.load() >= options.maxJobs);
    };
    const auto record = [&](const QueueRunResult& result, bool failed) {
        std::lock_guard<std::mutex> guard(lock);
        lastWork = std::chrono::steady_clock::now();
        if (failed) { ++summary.failed; return; }
        switch (result.state) {
        case QueueRunResult::State::completed:
            if (result.job && result.job->status == JobStatus::failed) ++summary.failed;
            else ++summary.completed;
            break;
        case QueueRunResult::State::cancelled: ++summary.cancelled; break;
        case QueueRunResult::State::requeued: ++summary.requeued; break;
        case QueueRunResult::State::paused: ++summary.paused; break;
        case QueueRunResult::State::empty: break;
        }
    };
    const auto worker = [&](std::size_t index) {
        const std::string name = options.worker + "-" + std::to_string(index);
        while (!stopping()) {
            if (index == 0) {
                const auto findings = watchdog(options.watchdogGraceSeconds);
                std::lock_guard<std::mutex> guard(lock);
                for (const auto& finding : findings) summary.quarantined += finding.quarantined;
            }
            // Reserve a slot before claiming so maxJobs is never exceeded by racing workers.
            if (options.maxJobs && started.fetch_add(1) >= options.maxJobs) break;
            if (!options.maxJobs) started.fetch_add(1);
            QueueRunResult result;
            bool failed = false;
            try {
                if (options.commandTemplate.empty()) {
                    result = runNext(name, {}, options.limits, {}, options.admission);
                } else {
                    dynsyslab::core::JobStore store(database_);
                    const auto claimed = store.claimNext(name, options.admission,
                                                         dynsyslab::core::currentProcessId());
                    if (claimed) result = runRemote(*claimed, options.commandTemplate);
                }
            } catch (const std::exception&) { failed = true; }  // recorded on the job
            if (!failed && result.state == QueueRunResult::State::empty) {
                started.fetch_sub(1);
                bool idle;
                { std::lock_guard<std::mutex> guard(lock);
                  idle = options.idleExitSeconds >= 0.0
                      && std::chrono::duration<double>(std::chrono::steady_clock::now() - lastWork)
                             .count() >= options.idleExitSeconds; }
                if (idle) { finished = true; break; }
                std::this_thread::sleep_for(std::chrono::duration<double>(options.pollSeconds));
                continue;
            }
            record(result, failed);
        }
    };
    std::vector<std::thread> threads;
    for (std::size_t index = 1; index < options.workers; ++index) threads.emplace_back(worker, index);
    worker(0);
    finished = true;
    for (auto& thread : threads) thread.join();
    return summary;
}

}  // namespace dynsyslab::workflow
