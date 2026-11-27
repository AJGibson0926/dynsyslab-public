// Claim-time scheduling for the native job store: dependency conditions, job-array caps,
// spec claims and admission (Slurm semantics; Yoo, Jette & Grondona 2003), plus the pause
// watchdog that quarantines workers which never acknowledge a pause.
#include "JobStoreInternal.hpp"

#include <algorithm>
#include <cerrno>
#include <cstdlib>
#include <map>
#include <stdexcept>

#ifdef _WIN32
#include <windows.h>
#else
#include <csignal>
#include <unistd.h>
#endif

namespace dynsyslab::core {

using namespace jobstore;

std::string dependencyConditionName(DependencyCondition condition) {
    switch (condition) {
    case DependencyCondition::afterok: return "afterok";
    case DependencyCondition::afterany: return "afterany";
    case DependencyCondition::afternotok: return "afternotok";
    case DependencyCondition::aftercorr: return "aftercorr";
    }
    return "afterok";
}

DependencyCondition parseDependencyCondition(const std::string& value) {
    for (const auto condition : {DependencyCondition::afterok, DependencyCondition::afterany,
                                 DependencyCondition::afternotok, DependencyCondition::aftercorr})
        if (dependencyConditionName(condition) == value) return condition;
    throw std::invalid_argument("dependency condition must be afterok, afterany, afternotok or aftercorr");
}

DependencyState dependencyState(DependencyCondition condition, std::optional<JobStatus> status) {
    if (!status) return DependencyState::impossible;  // the named job does not exist
    const bool ok = *status == JobStatus::completed;
    const bool notOk = *status == JobStatus::failed || *status == JobStatus::cancelled
        || *status == JobStatus::quarantined;
    switch (condition) {
    case DependencyCondition::afterok:
    case DependencyCondition::aftercorr:
        return ok ? DependencyState::satisfied
                  : notOk ? DependencyState::impossible : DependencyState::waiting;
    case DependencyCondition::afterany:
        return ok || notOk ? DependencyState::satisfied : DependencyState::waiting;
    case DependencyCondition::afternotok:
        return notOk ? DependencyState::satisfied
                     : ok ? DependencyState::impossible : DependencyState::waiting;
    }
    return DependencyState::waiting;
}

std::string admissionDenial(const JobResources& candidate, const std::vector<JobRecord>& running,
                            const AdmissionPolicy& policy) {
    if (policy.maxRunning && running.size() >= policy.maxRunning) return "max_running";
    std::size_t threads = std::max<std::size_t>(1, candidate.threads);
    std::size_t memory = candidate.memoryBytes;
    for (const auto& job : running) {
        threads += std::max<std::size_t>(1, job.resources.threads);
        memory += job.resources.memoryBytes;
    }
    // A single job larger than the whole pool is still admitted when nothing else runs, so a
    // limit never deadlocks the queue; the job's own limits enforce it during execution.
    if (running.empty()) return {};
    if (policy.maxThreads && threads > policy.maxThreads) return "max_threads";
    if (policy.maxMemoryBytes && memory > policy.maxMemoryBytes) return "max_memory";
    return {};
}

namespace {

std::optional<JobStatus> statusOf(const std::map<std::string, JobRecord>& byId,
                                  const std::string& id) {
    const auto found = byId.find(id);
    if (found == byId.end()) return std::nullopt;
    return found->second.status;
}

// Evaluates all dependencies of one job; the first impossible one is returned in `reason`.
DependencyState dependenciesOf(const JobRecord& job, const std::map<std::string, JobRecord>& byId,
                               std::string& reason) {
    DependencyState overall = DependencyState::satisfied;
    for (const auto& dependency : job.dependencies) {
        std::optional<JobStatus> status;
        if (dependency.condition == DependencyCondition::aftercorr) {
            for (const auto& item : byId)
                if (item.second.arrayId == dependency.job
                    && item.second.arrayIndex == job.arrayIndex) status = item.second.status;
        } else status = statusOf(byId, dependency.job);
        const auto state = dependencyState(dependency.condition, status);
        if (state == DependencyState::impossible) {
            reason = "dependency " + dependencyConditionName(dependency.condition) + ':'
                + dependency.job + (status ? " ended " + jobStatusName(*status) : " does not exist");
            return state;
        }
        if (state == DependencyState::waiting) overall = DependencyState::waiting;
    }
    return overall;
}

}  // namespace

std::optional<JobRecord> JobStore::claimNext(const std::string& worker,
                                             const AdmissionPolicy& policy,
                                             std::uint64_t processId) {
    auto* database = implementation_->database;
    execute(database, "BEGIN IMMEDIATE");
    try {
        const auto queued = select(database, "WHERE status='queued' ORDER BY rowid");
        if (queued.empty()) { execute(database, "COMMIT"); return std::nullopt; }
        const auto all = select(database, "");
        std::map<std::string, JobRecord> byId;
        std::vector<JobRecord> running;
        for (const auto& job : all) {
            byId.emplace(job.id, job);
            if (job.status == JobStatus::running) running.push_back(job);
        }
        std::optional<JobRecord> chosen;
        for (const auto& job : queued) {
            std::string reason;
            const auto dependencies = dependenciesOf(job, byId, reason);
            if (dependencies == DependencyState::impossible) {
                setStatus(database, job.id, JobStatus::cancelled, reason, 0.0);
                continue;
            }
            if (dependencies == DependencyState::waiting) continue;
            if (!job.arrayId.empty() && job.arrayMaxConcurrent
                && static_cast<std::size_t>(std::count_if(running.begin(), running.end(),
                       [&](const JobRecord& other) { return other.arrayId == job.arrayId; }))
                       >= job.arrayMaxConcurrent) continue;
            if (!job.specKey.empty()
                && std::any_of(running.begin(), running.end(), [&](const JobRecord& other) {
                       return other.specKey == job.specKey; })) continue;
            if (!admissionDenial(job.resources, running, policy).empty()) continue;
            chosen = job;
            break;
        }
        if (chosen) {
            const std::string now = timestamp();
            Statement(database,
                "UPDATE jobs SET status='running',worker=?,attempts=attempts+1,updated_utc=?,"
                "process_id=?,host=?,pause_requested=0 WHERE id=? AND status='queued'", "claim job")
                .text(worker).text(now).integer(static_cast<long long>(processId))
                .text(processId ? currentHost() : std::string{}).text(chosen->id).run();
            chosen->status = JobStatus::running; chosen->worker = worker;
            chosen->updatedUtc = now; ++chosen->attempts; chosen->processId = processId;
        }
        execute(database, "COMMIT");
        return chosen;
    } catch (...) { execute(database, "ROLLBACK"); throw; }
}

std::vector<WatchdogFinding> JobStore::watchPausedJobs(double graceSeconds) {
    auto* database = implementation_->database;
    std::vector<WatchdogFinding> findings;
    Statement query(database,
        std::string("SELECT ") + columns + ",(julianday('now')-julianday(updated_utc))*86400.0"
        " FROM jobs WHERE status='running' AND pause_requested=1", "query paused jobs");
    std::vector<std::pair<JobRecord, double>> rows;
    while (query.row())
        rows.emplace_back(record(query.get()), sqlite3_column_double(query.get(), 22));
    const auto host = currentHost();
    for (const auto& [job, age] : rows) {
        WatchdogFinding finding;
        finding.job = job.id; finding.heartbeatAgeSeconds = age;
        // A worker on another host cannot be probed here; only its heartbeat is judged.
        const bool local = job.processId && job.host == host;
        finding.workerAlive = local ? processAlive(job.processId) : age <= graceSeconds;
        if (local && !finding.workerAlive) finding.classification = "exited_without_acknowledgement";
        else if (age > graceSeconds) finding.classification = "rogue_suspected";
        else finding.classification = "finishing_current_work";
        if (finding.classification != "finishing_current_work") {
            // Quarantine, not requeue: a live rogue worker must not be duplicated.
            setStatus(database, job.id, JobStatus::quarantined,
                      "pause watchdog: " + finding.classification, 0.0);
            finding.quarantined = true;
        }
        findings.push_back(std::move(finding));
    }
    return findings;
}

bool processAlive(std::uint64_t processId) {
    if (!processId) return false;
#ifdef _WIN32
    HANDLE handle = OpenProcess(PROCESS_QUERY_LIMITED_INFORMATION, FALSE,
                                static_cast<DWORD>(processId));
    if (!handle) return false;
    DWORD code = 0;
    const bool alive = GetExitCodeProcess(handle, &code) && code == STILL_ACTIVE;
    CloseHandle(handle);
    return alive;
#else
    return ::kill(static_cast<pid_t>(processId), 0) == 0 || errno == EPERM;
#endif
}

std::uint64_t currentProcessId() {
#ifdef _WIN32
    return GetCurrentProcessId();
#else
    return static_cast<std::uint64_t>(::getpid());
#endif
}

std::string currentHost() {
#ifdef _WIN32
    wchar_t name[MAX_COMPUTERNAME_LENGTH + 1]{};
    DWORD size = MAX_COMPUTERNAME_LENGTH + 1;
    if (!GetComputerNameW(name, &size)) return "localhost";
    std::string result;
    for (DWORD index = 0; index < size; ++index) result += static_cast<char>(name[index] & 0x7f);
    return result;
#else
    char name[256]{};
    return gethostname(name, sizeof(name)) == 0 ? std::string(name) : std::string("localhost");
#endif
}

}  // namespace dynsyslab::core
