#pragma once
// Shared SQLite helpers for the JobStore translation units.

#include <dynsyslab/core/JobStore.hpp>

#include <sqlite3.h>

#include <string>
#include <vector>

namespace dynsyslab::core {

struct JobStore::Implementation {
    sqlite3* database{};
    explicit Implementation(const std::filesystem::path& path);
    ~Implementation();
};

namespace jobstore {

void check(int status, sqlite3* database, const char* operation);
void execute(sqlite3* database, const char* sql);
std::string timestamp();
std::string column(sqlite3_stmt* statement, int index);

// One prepared statement, finalized on scope exit; text bindings are copied.
class Statement {
public:
    Statement(sqlite3* database, const std::string& sql, const char* operation);
    ~Statement();
    Statement(const Statement&) = delete;
    Statement& operator=(const Statement&) = delete;
    Statement& text(const std::string& value);
    Statement& integer(long long value);
    Statement& real(double value);
    bool row();   // sqlite3_step; true while a row is available
    void run();   // step to completion
    sqlite3_stmt* get() const { return statement_; }

private:
    sqlite3* database_{};
    sqlite3_stmt* statement_{};
    const char* operation_{};
    int parameter_{1};
};

extern const char* const columns;
JobRecord record(sqlite3_stmt* statement);
std::string encodeDependencies(const std::vector<JobDependency>& dependencies);
std::vector<JobDependency> decodeDependencies(const std::string& text);
// Sets status, message and progress unconditionally (terminal and administrative moves).
void setStatus(sqlite3* database, const std::string& id, JobStatus status,
               const std::string& message, double progress);
std::vector<JobRecord> select(sqlite3* database, const std::string& where,
                              const std::vector<std::string>& texts = {});

}  // namespace jobstore
}  // namespace dynsyslab::core
