#include <dynsyslab/core/JobStore.hpp>

#include <sqlite3.h>

#include <algorithm>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

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

std::string identifier(const std::string& manifest) {
    const auto ticks = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::uint64_t hash = 1469598103934665603ull ^ static_cast<std::uint64_t>(ticks);
    for (const unsigned char value : manifest) { hash ^= value; hash *= 1099511628211ull; }
    std::ostringstream out; out << std::hex << std::setw(16) << std::setfill('0') << hash;
    return out.str();
}

JobStatus parseStatus(const std::string& value) {
    if (value == "queued") return JobStatus::queued;
    if (value == "running") return JobStatus::running;
    if (value == "completed") return JobStatus::completed;
    if (value == "failed") return JobStatus::failed;
    if (value == "cancelled") return JobStatus::cancelled;
    throw std::runtime_error("unknown job status: " + value);
}

std::string column(sqlite3_stmt* statement, int index) {
    const auto* value = sqlite3_column_text(statement, index);
    return value ? reinterpret_cast<const char*>(value) : std::string{};
}

JobRecord record(sqlite3_stmt* statement) {
    JobRecord result;
    result.id = column(statement, 0); result.label = column(statement, 1);
    result.kind = column(statement, 2); result.manifestJson = column(statement, 3);
    result.status = parseStatus(column(statement, 4));
    result.progress = sqlite3_column_double(statement, 5); result.message = column(statement, 6);
    result.worker = column(statement, 7); result.createdUtc = column(statement, 8);
    result.updatedUtc = column(statement, 9);
    result.attempts = static_cast<std::size_t>(sqlite3_column_int64(statement, 10));
    return result;
}

constexpr const char* columns =
    "id,label,kind,manifest_json,status,progress,message,worker,created_utc,updated_utc,attempts";

}  // namespace

struct JobStore::Implementation {
    sqlite3* database{};
    explicit Implementation(const std::filesystem::path& path) {
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
    }
    ~Implementation() { if (database) sqlite3_close(database); }
};

JobStore::JobStore(const std::filesystem::path& path)
    : implementation_(std::make_unique<Implementation>(path)) {}
JobStore::~JobStore() = default;

std::string jobStatusName(JobStatus status) {
    switch (status) {
    case JobStatus::queued: return "queued"; case JobStatus::running: return "running";
    case JobStatus::completed: return "completed"; case JobStatus::failed: return "failed";
    case JobStatus::cancelled: return "cancelled";
    }
    return "unknown";
}

std::string JobStore::submit(const std::string& label, const std::string& kind,
                             const std::string& manifestJson) {
    const std::string id = identifier(manifestJson), now = timestamp();
    sqlite3_stmt* statement = nullptr;
    check(sqlite3_prepare_v2(implementation_->database,
        "INSERT INTO jobs(id,label,kind,manifest_json,status,created_utc,updated_utc)"
        " VALUES(?,?,?,?,'queued',?,?)", -1, &statement, nullptr),
        implementation_->database, "prepare job insert");
    sqlite3_bind_text(statement, 1, id.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(statement, 2, label.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(statement, 3, kind.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(statement, 4, manifestJson.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(statement, 5, now.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(statement, 6, now.c_str(), -1, SQLITE_TRANSIENT);
    const int status = sqlite3_step(statement); sqlite3_finalize(statement);
    check(status, implementation_->database, "insert job");
    return id;
}

std::optional<JobRecord> JobStore::claimNext(const std::string& worker) {
    execute(implementation_->database, "BEGIN IMMEDIATE");
    try {
        sqlite3_stmt* query = nullptr;
        const std::string sql = std::string("SELECT ") + columns
            + " FROM jobs WHERE status='queued' ORDER BY created_utc,id LIMIT 1";
        check(sqlite3_prepare_v2(implementation_->database, sql.c_str(), -1, &query, nullptr),
              implementation_->database, "prepare queued job query");
        if (sqlite3_step(query) != SQLITE_ROW) {
            sqlite3_finalize(query); execute(implementation_->database, "COMMIT"); return std::nullopt;
        }
        JobRecord result = record(query); sqlite3_finalize(query);
        sqlite3_stmt* update = nullptr;
        check(sqlite3_prepare_v2(implementation_->database,
            "UPDATE jobs SET status='running',worker=?,attempts=attempts+1,updated_utc=?"
            " WHERE id=? AND status='queued'", -1, &update, nullptr),
            implementation_->database, "prepare job claim");
        const std::string now = timestamp();
        sqlite3_bind_text(update, 1, worker.c_str(), -1, SQLITE_TRANSIENT);
        sqlite3_bind_text(update, 2, now.c_str(), -1, SQLITE_TRANSIENT);
        sqlite3_bind_text(update, 3, result.id.c_str(), -1, SQLITE_TRANSIENT);
        check(sqlite3_step(update), implementation_->database, "claim job");
        sqlite3_finalize(update); execute(implementation_->database, "COMMIT");
        result.status = JobStatus::running; result.worker = worker;
        result.updatedUtc = now; ++result.attempts; return result;
    } catch (...) { execute(implementation_->database, "ROLLBACK"); throw; }
}

void JobStore::updateProgress(const std::string& id, double progress,
                              const std::string& message) {
    sqlite3_stmt* statement = nullptr;
    check(sqlite3_prepare_v2(implementation_->database,
        "UPDATE jobs SET progress=?,message=?,updated_utc=? WHERE id=? AND status='running'",
        -1, &statement, nullptr), implementation_->database, "prepare progress update");
    const std::string now = timestamp();
    sqlite3_bind_double(statement, 1, std::clamp(progress, 0.0, 1.0));
    sqlite3_bind_text(statement, 2, message.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(statement, 3, now.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(statement, 4, id.c_str(), -1, SQLITE_TRANSIENT);
    check(sqlite3_step(statement), implementation_->database, "update progress");
    sqlite3_finalize(statement);
}

void terminal(sqlite3* database, const std::string& id, JobStatus status,
              const std::string& message) {
    sqlite3_stmt* statement = nullptr;
    check(sqlite3_prepare_v2(database,
        "UPDATE jobs SET status=?,progress=?,message=?,updated_utc=? WHERE id=?",
        -1, &statement, nullptr), database, "prepare terminal job update");
    const std::string now = timestamp(), name = jobStatusName(status);
    sqlite3_bind_text(statement, 1, name.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_double(statement, 2, status == JobStatus::completed ? 1.0 : 0.0);
    sqlite3_bind_text(statement, 3, message.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(statement, 4, now.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(statement, 5, id.c_str(), -1, SQLITE_TRANSIENT);
    check(sqlite3_step(statement), database, "update terminal job"); sqlite3_finalize(statement);
}

void JobStore::complete(const std::string& id, const std::string& message) {
    terminal(implementation_->database, id, JobStatus::completed, message);
}
void JobStore::fail(const std::string& id, const std::string& message) {
    terminal(implementation_->database, id, JobStatus::failed, message);
}
void JobStore::cancel(const std::string& id, const std::string& message) {
    terminal(implementation_->database, id, JobStatus::cancelled, message);
}

std::size_t JobStore::recoverStale(double ageSeconds) {
    sqlite3_stmt* statement = nullptr;
    check(sqlite3_prepare_v2(implementation_->database,
        "UPDATE jobs SET status='queued',worker='',message='recovered stale worker',updated_utc=?"
        " WHERE status='running' AND (julianday('now')-julianday(updated_utc))*86400.0>?",
        -1, &statement, nullptr), implementation_->database, "prepare stale recovery");
    const std::string now = timestamp();
    sqlite3_bind_text(statement, 1, now.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_double(statement, 2, std::max(0.0, ageSeconds));
    check(sqlite3_step(statement), implementation_->database, "recover stale jobs");
    const auto count = static_cast<std::size_t>(sqlite3_changes(implementation_->database));
    sqlite3_finalize(statement); return count;
}

std::vector<JobRecord> JobStore::list(std::size_t limit,
                                      std::optional<JobStatus> status) const {
    sqlite3_stmt* statement = nullptr;
    std::string sql = std::string("SELECT ") + columns + " FROM jobs";
    if (status) sql += " WHERE status=?";
    sql += " ORDER BY created_utc DESC,id DESC LIMIT ?";
    check(sqlite3_prepare_v2(implementation_->database, sql.c_str(), -1, &statement, nullptr),
          implementation_->database, "prepare job list");
    int parameter = 1;
    std::string name;
    if (status) { name = jobStatusName(*status);
        sqlite3_bind_text(statement, parameter++, name.c_str(), -1, SQLITE_TRANSIENT); }
    sqlite3_bind_int64(statement, parameter, static_cast<sqlite3_int64>(limit));
    std::vector<JobRecord> result;
    while (sqlite3_step(statement) == SQLITE_ROW) result.push_back(record(statement));
    sqlite3_finalize(statement); return result;
}

std::optional<JobRecord> JobStore::find(const std::string& id) const {
    sqlite3_stmt* statement = nullptr;
    const std::string sql = std::string("SELECT ") + columns + " FROM jobs WHERE id=?";
    check(sqlite3_prepare_v2(implementation_->database, sql.c_str(), -1, &statement, nullptr),
          implementation_->database, "prepare job lookup");
    sqlite3_bind_text(statement, 1, id.c_str(), -1, SQLITE_TRANSIENT);
    if (sqlite3_step(statement) != SQLITE_ROW) { sqlite3_finalize(statement); return std::nullopt; }
    JobRecord result = record(statement); sqlite3_finalize(statement); return result;
}

}  // namespace dynsyslab::core
