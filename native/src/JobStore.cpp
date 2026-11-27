#include "JobStoreInternal.hpp"

#include <algorithm>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <set>
#include <sstream>
#include <stdexcept>

namespace dynsyslab::core {
namespace jobstore {

void check(int status, sqlite3* database, const char* operation) {
    if (status != SQLITE_OK && status != SQLITE_DONE && status != SQLITE_ROW)
        throw std::runtime_error(std::string(operation) + ": " + sqlite3_errmsg(database));
}

void execute(sqlite3* database, const char* sql) {
    char* message = nullptr;
    const int status = sqlite3_exec(database, sql, nullptr, nullptr, &message);
    if (status != SQLITE_OK) {
        const std::string text = message ? message : sqlite3_errmsg(database);
        sqlite3_free(message);
        throw std::runtime_error("SQLite: " + text);
    }
}

std::string timestamp() {
    const auto now = std::chrono::system_clock::now();
    const std::time_t raw = std::chrono::system_clock::to_time_t(now);
    std::tm utc{};
#ifdef _WIN32
    gmtime_s(&utc, &raw);
#else
    gmtime_r(&raw, &utc);
#endif
    std::ostringstream out; out << std::put_time(&utc, "%Y-%m-%dT%H:%M:%SZ");
    return out.str();
}

std::string column(sqlite3_stmt* statement, int index) {
    const auto* value = sqlite3_column_text(statement, index);
    return value ? reinterpret_cast<const char*>(value) : std::string{};
}

Statement::Statement(sqlite3* database, const std::string& sql, const char* operation)
    : database_(database), operation_(operation) {
    check(sqlite3_prepare_v2(database, sql.c_str(), -1, &statement_, nullptr), database, operation);
}
Statement::~Statement() { sqlite3_finalize(statement_); }
Statement& Statement::text(const std::string& value) {
    sqlite3_bind_text(statement_, parameter_++, value.c_str(), -1, SQLITE_TRANSIENT); return *this;
}
Statement& Statement::integer(long long value) {
    sqlite3_bind_int64(statement_, parameter_++, static_cast<sqlite3_int64>(value)); return *this;
}
Statement& Statement::real(double value) {
    sqlite3_bind_double(statement_, parameter_++, value); return *this;
}
bool Statement::row() {
    const int status = sqlite3_step(statement_);
    check(status, database_, operation_);
    return status == SQLITE_ROW;
}
void Statement::run() { while (row()) {} }

const char* const columns =
    "id,label,kind,manifest_json,status,progress,message,worker,created_utc,updated_utc,attempts,"
    "retry_limit,dependencies,threads,memory_bytes,array_id,array_index,array_max_concurrent,"
    "spec_key,process_id,host,pause_requested";

JobRecord record(sqlite3_stmt* statement) {
    const auto count = [statement](int index) {
        return static_cast<std::size_t>(sqlite3_column_int64(statement, index));
    };
    JobRecord result;
    result.id = column(statement, 0); result.label = column(statement, 1);
    result.kind = column(statement, 2); result.manifestJson = column(statement, 3);
    result.status = parseJobStatus(column(statement, 4));
    result.progress = sqlite3_column_double(statement, 5); result.message = column(statement, 6);
    result.worker = column(statement, 7); result.createdUtc = column(statement, 8);
    result.updatedUtc = column(statement, 9); result.attempts = count(10);
    result.retryLimit = count(11); result.dependencies = decodeDependencies(column(statement, 12));
    result.resources = {count(13), count(14)};
    result.arrayId = column(statement, 15); result.arrayIndex = count(16);
    result.arrayMaxConcurrent = count(17); result.specKey = column(statement, 18);
    result.processId = static_cast<std::uint64_t>(sqlite3_column_int64(statement, 19));
    result.host = column(statement, 20); result.pauseRequested = count(21) != 0;
    return result;
}

std::string encodeDependencies(const std::vector<JobDependency>& dependencies) {
    std::string result;
    for (const auto& item : dependencies) {
        if (item.job.empty() || item.job.find_first_of(":;") != std::string::npos)
            throw std::invalid_argument("dependency job ids must be nonempty and contain no ':' or ';'");
        if (!result.empty()) result += ';';
        result += dependencyConditionName(item.condition) + ':' + item.job;
    }
    return result;
}

std::vector<JobDependency> decodeDependencies(const std::string& text) {
    std::vector<JobDependency> result;
    std::istringstream stream(text);
    for (std::string item; std::getline(stream, item, ';');) {
        const auto colon = item.find(':');
        if (colon == std::string::npos) continue;
        result.push_back({item.substr(colon + 1), parseDependencyCondition(item.substr(0, colon))});
    }
    return result;
}

void setStatus(sqlite3* database, const std::string& id, JobStatus status,
               const std::string& message, double progress) {
    Statement(database,
        "UPDATE jobs SET status=?,progress=?,message=?,updated_utc=?,pause_requested=0 WHERE id=?",
        "update job status")
        .text(jobStatusName(status)).real(progress).text(message).text(timestamp()).text(id).run();
}

std::vector<JobRecord> select(sqlite3* database, const std::string& where,
                              const std::vector<std::string>& texts) {
    Statement query(database, std::string("SELECT ") + columns + " FROM jobs " + where,
                    "query jobs");
    for (const auto& text : texts) query.text(text);
    std::vector<JobRecord> result;
    while (query.row()) result.push_back(record(query.get()));
    return result;
}

}  // namespace jobstore

using namespace jobstore;

namespace {

std::string identifier(const std::string& manifest) {
    const auto ticks = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::uint64_t hash = 1469598103934665603ull ^ static_cast<std::uint64_t>(ticks);
    for (const unsigned char value : manifest) { hash ^= value; hash *= 1099511628211ull; }
    std::ostringstream out; out << std::hex << std::setw(16) << std::setfill('0') << hash;
    return out.str();
}

// Columns added after the first native schema; older stores gain them in place.
void migrate(sqlite3* database) {
    std::set<std::string> present;
    {
        Statement info(database, "PRAGMA table_info(jobs)", "read job columns");
        while (info.row()) present.insert(column(info.get(), 1));
    }
    const char* added[][2] = {
        {"retry_limit", "INTEGER NOT NULL DEFAULT 0"}, {"dependencies", "TEXT NOT NULL DEFAULT ''"},
        {"threads", "INTEGER NOT NULL DEFAULT 0"}, {"memory_bytes", "INTEGER NOT NULL DEFAULT 0"},
        {"array_id", "TEXT NOT NULL DEFAULT ''"}, {"array_index", "INTEGER NOT NULL DEFAULT 0"},
        {"array_max_concurrent", "INTEGER NOT NULL DEFAULT 0"},
        {"spec_key", "TEXT NOT NULL DEFAULT ''"}, {"process_id", "INTEGER NOT NULL DEFAULT 0"},
        {"host", "TEXT NOT NULL DEFAULT ''"}, {"pause_requested", "INTEGER NOT NULL DEFAULT 0"}};
    for (const auto& item : added)
        if (!present.count(item[0]))
            execute(database, (std::string("ALTER TABLE jobs ADD COLUMN ") + item[0] + ' '
                               + item[1]).c_str());
    execute(database,
        "CREATE INDEX IF NOT EXISTS jobs_spec_key ON jobs(spec_key,status);"
        "CREATE INDEX IF NOT EXISTS jobs_array ON jobs(array_id,array_index);");
}

void requireStatus(sqlite3* database, const std::string& id,
                   std::initializer_list<JobStatus> allowed, const char* action) {
    const auto rows = select(database, "WHERE id=?", {id});
    if (rows.empty()) throw std::invalid_argument("queued job does not exist: " + id);
    for (const auto status : allowed) if (rows.front().status == status) return;
    throw std::invalid_argument(std::string(action) + " is not allowed for a "
                                + jobStatusName(rows.front().status) + " job");
}

}  // namespace

JobStore::Implementation::Implementation(const std::filesystem::path& path) {
    std::filesystem::create_directories(path.parent_path());
    check(sqlite3_open_v2(path.u8string().c_str(), &database,
        SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE | SQLITE_OPEN_FULLMUTEX, nullptr),
        database, "open job store");
    sqlite3_busy_timeout(database, 5000);
    execute(database, "PRAGMA journal_mode=WAL; PRAGMA synchronous=NORMAL;");
    execute(database,
        "CREATE TABLE IF NOT EXISTS jobs("
        "id TEXT PRIMARY KEY,label TEXT NOT NULL,kind TEXT NOT NULL,manifest_json TEXT NOT NULL,"
        "status TEXT NOT NULL,progress REAL NOT NULL DEFAULT 0,message TEXT NOT NULL DEFAULT '',"
        "worker TEXT NOT NULL DEFAULT '',created_utc TEXT NOT NULL,updated_utc TEXT NOT NULL,"
        "attempts INTEGER NOT NULL DEFAULT 0);"
        "CREATE INDEX IF NOT EXISTS jobs_status_created ON jobs(status,created_utc);");
    migrate(database);
}
JobStore::Implementation::~Implementation() { if (database) sqlite3_close(database); }

JobStore::JobStore(const std::filesystem::path& path)
    : implementation_(std::make_unique<Implementation>(path)) {}
JobStore::~JobStore() = default;

std::string jobStatusName(JobStatus status) {
    switch (status) {
    case JobStatus::queued: return "queued"; case JobStatus::running: return "running";
    case JobStatus::completed: return "completed"; case JobStatus::failed: return "failed";
    case JobStatus::cancelled: return "cancelled"; case JobStatus::paused: return "paused";
    case JobStatus::quarantined: return "quarantined"; case JobStatus::reserved: return "reserved";
    }
    return "unknown";
}

JobStatus parseJobStatus(const std::string& value) {
    for (const auto status : {JobStatus::queued, JobStatus::running, JobStatus::completed,
                              JobStatus::failed, JobStatus::cancelled, JobStatus::paused,
                              JobStatus::quarantined, JobStatus::reserved})
        if (jobStatusName(status) == value) return status;
    throw std::runtime_error("unknown job status: " + value);
}

std::string JobStore::submit(const std::string& label, const std::string& kind,
                             const std::string& manifestJson) {
    JobSubmission submission;
    submission.label = label; submission.kind = kind; submission.manifestJson = manifestJson;
    return submit(submission);
}

std::string JobStore::submit(const JobSubmission& item) {
    const std::string id = identifier(item.manifestJson + item.arrayId
                                      + std::to_string(item.arrayIndex));
    const std::string now = timestamp();
    Statement(implementation_->database,
        "INSERT INTO jobs(id,label,kind,manifest_json,status,created_utc,updated_utc,retry_limit,"
        "dependencies,threads,memory_bytes,array_id,array_index,array_max_concurrent,spec_key)"
        " VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", "insert job")
        .text(id).text(item.label).text(item.kind).text(item.manifestJson)
        .text(item.reserved ? "reserved" : "queued").text(now).text(now)
        .integer(static_cast<long long>(item.retryLimit)).text(encodeDependencies(item.dependencies))
        .integer(static_cast<long long>(item.resources.threads))
        .integer(static_cast<long long>(item.resources.memoryBytes))
        .text(item.arrayId).integer(static_cast<long long>(item.arrayIndex))
        .integer(static_cast<long long>(item.arrayMaxConcurrent)).text(item.specKey).run();
    return id;
}

void JobStore::updateProgress(const std::string& id, double progress,
                              const std::string& message) {
    Statement(implementation_->database,
        "UPDATE jobs SET progress=?,message=?,updated_utc=? WHERE id=? AND status='running'",
        "update progress")
        .real(std::clamp(progress, 0.0, 1.0)).text(message).text(timestamp()).text(id).run();
}

void JobStore::complete(const std::string& id, const std::string& message) {
    setStatus(implementation_->database, id, JobStatus::completed, message, 1.0);
}

bool JobStore::fail(const std::string& id, const std::string& message) {
    const auto current = find(id);
    if (current && current->attempts <= current->retryLimit) {
        Statement(implementation_->database,
            "UPDATE jobs SET status='queued',worker='',process_id=0,progress=0,message=?,"
            "updated_utc=? WHERE id=?", "requeue failed job")
            .text("retry " + std::to_string(current->attempts) + " of "
                  + std::to_string(current->retryLimit) + " after failure: " + message)
            .text(timestamp()).text(id).run();
        return true;
    }
    setStatus(implementation_->database, id, JobStatus::failed, message, 0.0);
    return false;
}

void JobStore::cancel(const std::string& id, const std::string& message) {
    setStatus(implementation_->database, id, JobStatus::cancelled, message, 0.0);
}

JobStatus JobStore::pause(const std::string& id) {
    auto* database = implementation_->database;
    requireStatus(database, id, {JobStatus::queued, JobStatus::reserved, JobStatus::running,
                                 JobStatus::paused}, "pause");
    const auto current = find(id);
    if (current->status == JobStatus::running) {
        Statement(database, "UPDATE jobs SET pause_requested=1,updated_utc=? WHERE id=?",
                  "request pause").text(timestamp()).text(id).run();
        return JobStatus::running;
    }
    if (current->status != JobStatus::paused)
        setStatus(database, id, JobStatus::paused, "paused before execution", 0.0);
    return JobStatus::paused;
}

void JobStore::acknowledgePause(const std::string& id, const std::string& message) {
    Statement(implementation_->database,
        "UPDATE jobs SET status='paused',worker='',process_id=0,pause_requested=0,progress=0,"
        "message=?,updated_utc=? WHERE id=?", "acknowledge pause")
        .text(message).text(timestamp()).text(id).run();
}

void JobStore::resume(const std::string& id) {
    requireStatus(implementation_->database, id, {JobStatus::paused, JobStatus::quarantined},
                  "resume");
    Statement(implementation_->database,
        "UPDATE jobs SET status='queued',worker='',process_id=0,pause_requested=0,"
        "message='resumed',updated_utc=? WHERE id=?", "resume job")
        .text(timestamp()).text(id).run();
}

void JobStore::release(const std::string& id) {
    requireStatus(implementation_->database, id, {JobStatus::reserved}, "release");
    setStatus(implementation_->database, id, JobStatus::queued, "released from the side tray", 0.0);
}

void JobStore::remove(const std::string& id) {
    requireStatus(implementation_->database, id, {JobStatus::reserved}, "remove");
    Statement(implementation_->database, "DELETE FROM jobs WHERE id=? AND status='reserved'",
              "remove reserved job").text(id).run();
}

void JobStore::quarantine(const std::string& id, const std::string& message) {
    setStatus(implementation_->database, id, JobStatus::quarantined, message, 0.0);
}

std::size_t JobStore::recoverStale(double ageSeconds) {
    Statement statement(implementation_->database,
        "UPDATE jobs SET status='queued',worker='',process_id=0,message='recovered stale worker',"
        "updated_utc=? WHERE status='running' AND pause_requested=0"
        " AND (julianday('now')-julianday(updated_utc))*86400.0>?", "recover stale jobs");
    statement.text(timestamp()).real(std::max(0.0, ageSeconds)).run();
    return static_cast<std::size_t>(sqlite3_changes(implementation_->database));
}

std::vector<JobRecord> JobStore::list(std::size_t limit,
                                      std::optional<JobStatus> status) const {
    const std::string order = " ORDER BY created_utc DESC,rowid DESC LIMIT "
        + std::to_string(limit);
    return status ? select(implementation_->database, "WHERE status=?" + order,
                           {jobStatusName(*status)})
                  : select(implementation_->database, order);
}

std::optional<JobRecord> JobStore::find(const std::string& id) const {
    auto rows = select(implementation_->database, "WHERE id=?", {id});
    if (rows.empty()) return std::nullopt;
    return std::move(rows.front());
}

std::optional<JobRecord> JobStore::findBySpec(const std::string& specKey,
                                              JobStatus status) const {
    if (specKey.empty()) return std::nullopt;
    auto rows = select(implementation_->database,
        "WHERE spec_key=? AND status=? ORDER BY rowid DESC LIMIT 1",
        {specKey, jobStatusName(status)});
    if (rows.empty()) return std::nullopt;
    return std::move(rows.front());
}

std::vector<JobRecord> JobStore::arrayTasks(const std::string& arrayId) const {
    return select(implementation_->database, "WHERE array_id=? ORDER BY array_index", {arrayId});
}

}  // namespace dynsyslab::core
