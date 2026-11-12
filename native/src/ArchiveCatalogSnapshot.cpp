#include <dynsyslab/core/ArchiveIndex.hpp>

#include <sqlite3.h>

#include <stdexcept>

namespace dynsyslab::core {
namespace {

struct Database {
    sqlite3* value{};
    explicit Database(const std::filesystem::path& path, bool writable) {
        if (writable) std::filesystem::create_directories(path.parent_path());
        const int flags = (writable ? SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE
                                    : SQLITE_OPEN_READONLY) | SQLITE_OPEN_FULLMUTEX;
        if (sqlite3_open_v2(path.u8string().c_str(), &value, flags, nullptr) != SQLITE_OK)
            throw std::runtime_error("unable to open archive index snapshot store");
        sqlite3_busy_timeout(value, 10000);
        if (writable) {
            execute("CREATE TABLE IF NOT EXISTS archive_index_state("
                "singleton INTEGER PRIMARY KEY CHECK(singleton=1),"
                "generation INTEGER NOT NULL)");
            execute("INSERT OR IGNORE INTO archive_index_state(singleton,generation) "
                "VALUES(1,0)");
            execute("CREATE TABLE IF NOT EXISTS archive_catalog_snapshots("
                "consumer TEXT PRIMARY KEY,generation INTEGER NOT NULL,"
                "schema_version INTEGER NOT NULL,payload_json TEXT NOT NULL)");
        }
    }
    ~Database() { if (value) sqlite3_close(value); }
    void execute(const char* sql) {
        char* message{};
        if (sqlite3_exec(value, sql, nullptr, nullptr, &message) != SQLITE_OK) {
            const std::string text = message ? message : sqlite3_errmsg(value);
            sqlite3_free(message); throw std::runtime_error(text);
        }
    }
};

bool hasTable(Database& database, const char* name) {
    sqlite3_stmt* statement{};
    if (sqlite3_prepare_v2(database.value,
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?", -1,
        &statement, nullptr) != SQLITE_OK) return false;
    sqlite3_bind_text(statement, 1, name, -1, SQLITE_STATIC);
    const bool found = sqlite3_step(statement) == SQLITE_ROW;
    sqlite3_finalize(statement); return found;
}

std::uint64_t generation(Database& database) {
    sqlite3_stmt* statement{};
    if (sqlite3_prepare_v2(database.value,
        "SELECT generation FROM archive_index_state WHERE singleton=1", -1,
        &statement, nullptr) != SQLITE_OK) return 0;
    const auto result = sqlite3_step(statement) == SQLITE_ROW
        ? static_cast<std::uint64_t>(sqlite3_column_int64(statement, 0)) : 0;
    sqlite3_finalize(statement); return result;
}

}  // namespace

std::optional<std::string> loadArchiveCatalogSnapshot(
    const std::filesystem::path& archiveRoot, const std::string& consumer,
    std::uint32_t schemaVersion) {
    if (!std::filesystem::is_regular_file(archiveIndexPath(archiveRoot)))
        return std::nullopt;
    Database database(archiveIndexPath(archiveRoot), false);
    if (!hasTable(database, "archive_index_state") ||
        !hasTable(database, "archive_catalog_snapshots")) return std::nullopt;
    sqlite3_stmt* statement{};
    const char* sql =
        "SELECT snapshot.payload_json FROM archive_catalog_snapshots snapshot "
        "JOIN archive_index_state state ON state.singleton=1 "
        "AND snapshot.generation=state.generation "
        "WHERE snapshot.consumer=? AND snapshot.schema_version=?";
    if (sqlite3_prepare_v2(database.value, sql, -1, &statement, nullptr) != SQLITE_OK)
        return std::nullopt;
    sqlite3_bind_text(statement, 1, consumer.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_int(statement, 2, static_cast<int>(schemaVersion));
    std::optional<std::string> result;
    if (sqlite3_step(statement) == SQLITE_ROW) {
        const auto* text = sqlite3_column_text(statement, 0);
        if (text) result = reinterpret_cast<const char*>(text);
    }
    sqlite3_finalize(statement); return result;
}

bool storeArchiveCatalogSnapshot(
    const std::filesystem::path& archiveRoot, const std::string& consumer,
    std::uint32_t schemaVersion, std::uint64_t sourceGeneration,
    const std::string& payloadJson) {
    Database database(archiveIndexPath(archiveRoot), true);
    database.execute("BEGIN IMMEDIATE");
    if (generation(database) != sourceGeneration) {
        database.execute("ROLLBACK"); return false;
    }
    sqlite3_stmt* statement{};
    const char* sql = "INSERT INTO archive_catalog_snapshots(consumer,generation,"
        "schema_version,payload_json) VALUES(?,?,?,?) ON CONFLICT(consumer) DO "
        "UPDATE SET generation=excluded.generation,schema_version=excluded.schema_version,"
        "payload_json=excluded.payload_json";
    if (sqlite3_prepare_v2(database.value, sql, -1, &statement, nullptr) != SQLITE_OK) {
        database.execute("ROLLBACK"); throw std::runtime_error(sqlite3_errmsg(database.value));
    }
    sqlite3_bind_text(statement, 1, consumer.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_int64(statement, 2, static_cast<sqlite3_int64>(sourceGeneration));
    sqlite3_bind_int(statement, 3, static_cast<int>(schemaVersion));
    sqlite3_bind_text(statement, 4, payloadJson.c_str(), -1, SQLITE_TRANSIENT);
    const int status = sqlite3_step(statement); sqlite3_finalize(statement);
    if (status != SQLITE_DONE) {
        database.execute("ROLLBACK"); throw std::runtime_error(sqlite3_errmsg(database.value));
    }
    database.execute("COMMIT"); return true;
}

std::uint64_t archiveIndexGeneration(const std::filesystem::path& archiveRoot) {
    if (!std::filesystem::is_regular_file(archiveIndexPath(archiveRoot))) return 0;
    Database database(archiveIndexPath(archiveRoot), false);
    return hasTable(database, "archive_index_state") ? generation(database) : 0;
}

}  // namespace dynsyslab::core
