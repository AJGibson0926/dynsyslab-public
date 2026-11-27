#include <dynsyslab/core/ArchiveIndex.hpp>

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/FileDigest.hpp>
#include <dynsyslab/core/Json.hpp>

#include <sqlite3.h>

#include <map>
#include <set>
#include <stdexcept>

#ifdef _WIN32
#include <windows.h>
#endif

namespace dynsyslab::core {
namespace {

struct Database {
    sqlite3* value{};
    explicit Database(const std::filesystem::path& path, bool writable = true) {
        if (writable) std::filesystem::create_directories(path.parent_path());
        const int flags = (writable ? SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE
                                    : SQLITE_OPEN_READONLY) | SQLITE_OPEN_FULLMUTEX;
        if (sqlite3_open_v2(path.u8string().c_str(), &value, flags, nullptr) != SQLITE_OK)
            throw std::runtime_error("unable to open archive index");
        sqlite3_busy_timeout(value, 10000);
        if (!writable) return;
        execute("PRAGMA journal_mode=WAL");
        execute("PRAGMA synchronous=NORMAL");
        execute("CREATE TABLE IF NOT EXISTS archive_objects("
            "relative_path TEXT PRIMARY KEY, modified_ticks INTEGER NOT NULL,"
            "size_bytes INTEGER NOT NULL, sha256 TEXT NOT NULL,"
            "arrays_json TEXT NOT NULL,"
            "summary_version INTEGER NOT NULL DEFAULT 0,"
            "summary_json TEXT NOT NULL DEFAULT '')");
        if (!hasColumn("archive_objects", "summary_version"))
            execute("ALTER TABLE archive_objects ADD COLUMN "
                    "summary_version INTEGER NOT NULL DEFAULT 0");
        if (!hasColumn("archive_objects", "summary_json"))
            execute("ALTER TABLE archive_objects ADD COLUMN "
                    "summary_json TEXT NOT NULL DEFAULT ''");
        if (hasColumn("archive_objects", "metadata_json")) compactLegacyTable();
        execute("CREATE TABLE IF NOT EXISTS archive_index_state("
            "singleton INTEGER PRIMARY KEY CHECK(singleton=1),"
            "generation INTEGER NOT NULL)");
        execute("INSERT OR IGNORE INTO archive_index_state(singleton,generation) "
            "VALUES(1,0)");
        execute("CREATE TABLE IF NOT EXISTS archive_catalog_snapshots("
            "consumer TEXT PRIMARY KEY,generation INTEGER NOT NULL,"
            "schema_version INTEGER NOT NULL,payload_json TEXT NOT NULL)");
    }
    ~Database() { if (value) sqlite3_close(value); }
    void execute(const char* sql) {
        char* message{};
        if (sqlite3_exec(value, sql, nullptr, nullptr, &message) != SQLITE_OK) {
            const std::string text = message ? message : sqlite3_errmsg(value);
            sqlite3_free(message); throw std::runtime_error(text);
        }
    }
    bool hasColumn(const char* table, const char* column) {
        sqlite3_stmt* statement{};
        const std::string sql = std::string("PRAGMA table_info(") + table + ')';
        if (sqlite3_prepare_v2(value, sql.c_str(), -1, &statement, nullptr) != SQLITE_OK)
            return false;
        bool found = false;
        while (sqlite3_step(statement) == SQLITE_ROW) {
            const auto* name = sqlite3_column_text(statement, 1);
            if (name && std::string(reinterpret_cast<const char*>(name)) == column) {
                found = true; break;
            }
        }
        sqlite3_finalize(statement); return found;
    }
    bool hasTable(const char* table) {
        sqlite3_stmt* statement{};
        if (sqlite3_prepare_v2(value,
            "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?", -1,
            &statement, nullptr) != SQLITE_OK) return false;
        sqlite3_bind_text(statement, 1, table, -1, SQLITE_STATIC);
        const bool found = sqlite3_step(statement) == SQLITE_ROW;
        sqlite3_finalize(statement); return found;
    }
    void compactLegacyTable() {
        execute("BEGIN IMMEDIATE");
        try {
            execute("DROP TABLE IF EXISTS archive_objects_compact");
            execute("CREATE TABLE archive_objects_compact("
                "relative_path TEXT PRIMARY KEY, modified_ticks INTEGER NOT NULL,"
                "size_bytes INTEGER NOT NULL, sha256 TEXT NOT NULL,"
                "arrays_json TEXT NOT NULL,summary_version INTEGER NOT NULL DEFAULT 0,"
                "summary_json TEXT NOT NULL DEFAULT '')");
            execute("INSERT INTO archive_objects_compact SELECT relative_path,"
                "modified_ticks,size_bytes,sha256,arrays_json,summary_version,"
                "summary_json FROM archive_objects");
            execute("DROP TABLE archive_objects");
            execute("ALTER TABLE archive_objects_compact RENAME TO archive_objects");
            execute("COMMIT");
        } catch (...) {
            execute("ROLLBACK"); throw;
        }
        try {
            execute("PRAGMA wal_checkpoint(TRUNCATE)");
            execute("VACUUM");
        } catch (...) {
            // Compaction is an optimization. The committed schema migration is
            // valid even when another read-only process temporarily blocks VACUUM.
        }
    }
};

void copyField(json::Value::Object& target, const json::Value& source,
               const char* key) {
    const auto& value = source.at(key);
    if (!value.isNull()) target.emplace(key, value);
}

void copyPath(json::Value::Object& target, const json::Value& source,
              const char* outer, const char* middle, const char* inner = nullptr) {
    const auto& value = inner ? source.at(outer).at(middle).at(inner)
                              : source.at(outer).at(middle);
    if (value.isNull()) return;
    auto& outerObject = target[outer].object();
    if (!inner) outerObject.emplace(middle, value);
    else outerObject[middle].object().emplace(inner, value);
}

std::uint64_t modifiedTicks(const std::filesystem::path& path) {
#ifdef _WIN32
    auto native = std::filesystem::absolute(path).wstring();
    if (native.rfind(L"\\\\?\\", 0) != 0) native = L"\\\\?\\" + native;
    const HANDLE handle = CreateFileW(native.c_str(), FILE_READ_ATTRIBUTES,
        FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE, nullptr,
        OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, nullptr);
    if (handle == INVALID_HANDLE_VALUE)
        throw std::runtime_error("unable to read archive file identity");
    BY_HANDLE_FILE_INFORMATION data{};
    const bool available = GetFileInformationByHandle(handle, &data) != FALSE;
    CloseHandle(handle);
    if (!available) throw std::runtime_error("unable to read archive file identity");
    const std::uint64_t windowsTicks =
        (static_cast<std::uint64_t>(data.ftLastWriteTime.dwHighDateTime) << 32)
        | data.ftLastWriteTime.dwLowDateTime;
    constexpr std::uint64_t epochOffset = 116444736000000000ull;
    return windowsTicks >= epochOffset ? (windowsTicks - epochOffset) * 100ull : 0ull;
#else
    return static_cast<std::uint64_t>(
        std::filesystem::last_write_time(path).time_since_epoch().count());
#endif
}

std::string relativePath(const std::filesystem::path& root,
                         const std::filesystem::path& path) {
    const auto absoluteRoot = std::filesystem::absolute(root).lexically_normal();
    const auto absolutePath = std::filesystem::absolute(path).lexically_normal();
    const auto relative = absolutePath.lexically_relative(absoluteRoot);
    if (relative.empty() || *relative.begin() == L"..") return {};
    return relative.lexically_normal().generic_u8string();
}

std::string arrayJson(const std::vector<std::string>& arrays) {
    json::Value::Array values;
    for (const auto& name : arrays) values.emplace_back(name);
    return json::dump(json::Value(std::move(values)));
}

std::filesystem::path ioPath(const std::filesystem::path& path) {
#ifdef _WIN32
    auto native = std::filesystem::absolute(path).wstring();
    if (native.rfind(L"\\\\?\\", 0) != 0) native = L"\\\\?\\" + native;
    return std::filesystem::path(native);
#else
    return path;
#endif
}

std::vector<std::string> arrays(const std::string& encoded) {
    std::vector<std::string> result;
    // Keep the parsed document alive: iterating parse(...).array() directly would
    // range over a reference into a destroyed temporary.
    const auto parsed = json::parse(encoded);
    for (const auto& value : parsed.array())
        result.push_back(value.string());
    return result;
}

ArchiveIndexEntry describe(const std::filesystem::path& root,
                           const std::filesystem::path& path) {
    NpzIndex index;
    try { index = readNpzIndex(ioPath(path)); }
    catch (...) {
        // The pointer and exact digest are still valuable for formats that a
        // current decoder cannot inspect. Viewer keeps conservative labeling.
    }
    ArchiveIndexEntry result;
    result.path = std::filesystem::absolute(path).lexically_normal();
    result.relativePath = relativePath(root, result.path);
    result.modifiedTicks = modifiedTicks(result.path);
    result.sizeBytes = std::filesystem::file_size(ioPath(result.path));
    result.sha256 = sha256File(ioPath(result.path));
    result.arrays = index.arrays;
    result.summaryVersion = kArchiveSummaryVersion;
    result.summaryJson = summarizeArchiveMetadata(index.metadataJson);
    return result;
}

void upsert(Database& database, const ArchiveIndexEntry& entry) {
    sqlite3_stmt* statement{};
    const char* sql = "INSERT INTO archive_objects("
        "relative_path,modified_ticks,size_bytes,sha256,arrays_json,"
        "summary_version,summary_json) VALUES(?,?,?,?,?,?,?) "
        "ON CONFLICT(relative_path) DO UPDATE SET "
        "modified_ticks=excluded.modified_ticks,size_bytes=excluded.size_bytes,"
        "sha256=excluded.sha256,arrays_json=excluded.arrays_json,"
        "summary_version=excluded.summary_version,summary_json=excluded.summary_json";
    if (sqlite3_prepare_v2(database.value, sql, -1, &statement, nullptr) != SQLITE_OK)
        throw std::runtime_error(sqlite3_errmsg(database.value));
    sqlite3_bind_text(statement, 1, entry.relativePath.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_int64(statement, 2, static_cast<sqlite3_int64>(entry.modifiedTicks));
    sqlite3_bind_int64(statement, 3, static_cast<sqlite3_int64>(entry.sizeBytes));
    sqlite3_bind_text(statement, 4, entry.sha256.c_str(), -1, SQLITE_TRANSIENT);
    const auto names = arrayJson(entry.arrays);
    sqlite3_bind_text(statement, 5, names.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_int(statement, 6, static_cast<int>(entry.summaryVersion));
    sqlite3_bind_text(statement, 7, entry.summaryJson.c_str(), -1, SQLITE_TRANSIENT);
    const int status = sqlite3_step(statement); sqlite3_finalize(statement);
    if (status != SQLITE_DONE) throw std::runtime_error(sqlite3_errmsg(database.value));
}

void invalidateCatalogSnapshots(Database& database) {
    database.execute("UPDATE archive_index_state SET generation=generation+1 "
                     "WHERE singleton=1");
    database.execute("DELETE FROM archive_catalog_snapshots");
}

struct Signature {
    std::uint64_t modifiedTicks{};
    std::uint64_t sizeBytes{};
    std::uint32_t summaryVersion{};
    bool summaryReady{};
};

std::map<std::string, Signature> signatures(
    Database& database) {
    std::map<std::string, Signature> result;
    sqlite3_stmt* statement{};
    sqlite3_prepare_v2(database.value,
        "SELECT relative_path,modified_ticks,size_bytes,summary_version,"
        "length(summary_json) FROM archive_objects",
        -1, &statement, nullptr);
    while (sqlite3_step(statement) == SQLITE_ROW) {
        const auto* path = sqlite3_column_text(statement, 0);
        result.emplace(path ? reinterpret_cast<const char*>(path) : "", Signature{
            static_cast<std::uint64_t>(sqlite3_column_int64(statement, 1)),
            static_cast<std::uint64_t>(sqlite3_column_int64(statement, 2)),
            static_cast<std::uint32_t>(sqlite3_column_int(statement, 3)),
            sqlite3_column_int(statement, 4) > 0});
    }
    sqlite3_finalize(statement); return result;
}

bool archiveFile(const std::filesystem::path& path) {
    const auto name = path.filename().wstring();
    return path.extension() == L".npz" && name.find(L".partial") == std::wstring::npos;
}

std::size_t removeStaleSolutions(Database& database,
                                 const std::filesystem::path& archiveRoot) {
    if (!database.hasTable("solutions")) return 0;
    std::vector<std::string> stale;
    sqlite3_stmt* query{};
    if (sqlite3_prepare_v2(database.value, "SELECT path FROM solutions", -1,
                          &query, nullptr) != SQLITE_OK) return 0;
    while (sqlite3_step(query) == SQLITE_ROW) {
        const auto* text = sqlite3_column_text(query, 0);
        if (!text) continue;
        const std::string stored(reinterpret_cast<const char*>(text));
        auto path = std::filesystem::u8path(stored);
        if (path.is_relative()) path = archiveRoot / path;
        if (!std::filesystem::is_regular_file(ioPath(path))) stale.push_back(stored);
    }
    sqlite3_finalize(query);
    sqlite3_stmt* remove{};
    sqlite3_prepare_v2(database.value, "DELETE FROM solutions WHERE path=?", -1,
                       &remove, nullptr);
    std::size_t removed = 0;
    for (const auto& path : stale) {
        sqlite3_reset(remove); sqlite3_clear_bindings(remove);
        sqlite3_bind_text(remove, 1, path.c_str(), -1, SQLITE_TRANSIENT);
        if (sqlite3_step(remove) == SQLITE_DONE) ++removed;
    }
    sqlite3_finalize(remove); return removed;
}

}  // namespace

std::filesystem::path archiveIndexPath(const std::filesystem::path& archiveRoot) {
    return archiveRoot / L"management" / L"archive_index.sqlite";
}

std::string summarizeArchiveMetadata(const std::string& metadataJson) {
    if (metadataJson.empty()) return "{}";
    try {
        const auto source = json::parse(metadataJson);
        json::Value::Object summary;
        constexpr const char* fields[]{
            "label", "system_label", "system", "method", "analysis", "schema",
            "system_preset", "equations_version", "equations_fingerprint",
            "integration_form", "krom_mode_count", "eigenfunction_representation",
            "control_map_time_domain", "control_coordinate_system", "control_maps",
            "selected_modes", "independent_archive", "delta_t", "dataset_role",
            "attempt_count", "stop_reason", "controller_kind", "source_archive",
            "source_archives", "source_spectrum", "source_krom", "result_archives",
            "forcing_type", "data_role", "initial_state_policy", "dataset_references",
            "training_archives", "validation_archives", "test_archives", "validated_krom",
            "value_space", "eigenfunction"};
        for (const auto* field : fields) copyField(summary, source, field);
        copyPath(summary, source, "system_configuration", "identifier");
        copyPath(summary, source, "data_object", "kind");
        if (!source.at("analysis_metadata").at("krom").isNull())
            summary["analysis_metadata"].object().emplace("krom", true);
        copyPath(summary, source, "controller_configuration", "settings", "krom");
        copyPath(summary, source, "training_preflight", "source_krom");
        copyPath(summary, source, "controller", "controller");
        copyPath(summary, source, "controller", "control_representation");
        copyPath(summary, source, "controller", "target");
        copyPath(summary, source, "controller", "target_space");
        copyPath(summary, source, "controller", "rom", "source_archive");
        copyPath(summary, source, "system_parameters", "preset", "key");
        copyPath(summary, source, "system_parameters", "preset", "description");
        copyPath(summary, source, "system_parameters", "preset", "gas");
        copyPath(summary, source, "system_parameters", "preset", "temperature");
        copyPath(summary, source, "system_parameters", "preset", "temperature_c");
        return json::dump(json::Value(std::move(summary)));
    } catch (...) {
        return "{}";
    }
}

std::vector<ArchiveIndexEntry> loadArchiveIndex(
    const std::filesystem::path& archiveRoot) {
    if (!std::filesystem::is_regular_file(archiveIndexPath(archiveRoot))) return {};
    Database database(archiveIndexPath(archiveRoot), false);
    std::vector<ArchiveIndexEntry> result;
    sqlite3_stmt* statement{};
    const std::string sql = "SELECT relative_path,modified_ticks,size_bytes,sha256,"
        "arrays_json,summary_version,summary_json FROM archive_objects ORDER BY relative_path";
    if (sqlite3_prepare_v2(database.value, sql.c_str(), -1, &statement, nullptr) != SQLITE_OK)
        return result;
    while (sqlite3_step(statement) == SQLITE_ROW) {
        auto text = [&](int column) {
            const auto* value = sqlite3_column_text(statement, column);
            return std::string(value ? reinterpret_cast<const char*>(value) : "");
        };
        ArchiveIndexEntry entry;
        entry.relativePath = text(0);
        entry.path = archiveRoot / std::filesystem::u8path(entry.relativePath);
        entry.modifiedTicks = static_cast<std::uint64_t>(sqlite3_column_int64(statement, 1));
        entry.sizeBytes = static_cast<std::uint64_t>(sqlite3_column_int64(statement, 2));
        entry.sha256 = text(3);
        try { entry.arrays = arrays(text(4)); } catch (...) {}
        entry.summaryVersion = static_cast<std::uint32_t>(
            sqlite3_column_int(statement, 5));
        entry.summaryJson = text(6);
        result.push_back(std::move(entry));
    }
    sqlite3_finalize(statement); return result;
}

ArchiveIndexUpdate indexArchiveFiles(
    const std::filesystem::path& archiveRoot,
    const std::vector<std::filesystem::path>& paths) {
    ArchiveIndexUpdate result;
    if (archiveRoot.empty()) return result;
    Database database(archiveIndexPath(archiveRoot));
    database.execute("BEGIN IMMEDIATE");
    try {
        const auto known = signatures(database);
        for (auto path : paths) {
            if (path.is_relative()) path = archiveRoot / path;
            if (!std::filesystem::is_regular_file(ioPath(path)) || !archiveFile(path)) continue;
            ++result.discovered;
            const auto relative = relativePath(archiveRoot, path);
            if (relative.empty()) continue;
            const auto stamp = modifiedTicks(path);
            const auto bytes = std::filesystem::file_size(ioPath(path));
            const auto found = known.find(relative);
            if (found != known.end() && found->second.modifiedTicks == stamp &&
                found->second.sizeBytes == bytes && found->second.summaryReady &&
                found->second.summaryVersion == kArchiveSummaryVersion) {
                ++result.unchanged; continue;
            }
            upsert(database, describe(archiveRoot, path)); ++result.indexed;
        }
        if (result.indexed) invalidateCatalogSnapshots(database);
        database.execute("COMMIT");
    } catch (...) { database.execute("ROLLBACK"); throw; }
    return result;
}

ArchiveIndexUpdate reconcileArchiveIndex(
    const std::filesystem::path& archiveRoot) {
    ArchiveIndexUpdate result;
    if (!std::filesystem::is_directory(archiveRoot)) return result;
    Database database(archiveIndexPath(archiveRoot));
    auto known = signatures(database);
    std::set<std::string> present;
    database.execute("BEGIN IMMEDIATE");
    try {
        for (const auto& item : std::filesystem::recursive_directory_iterator(
                 archiveRoot, std::filesystem::directory_options::skip_permission_denied)) {
            if (!item.is_regular_file() || !archiveFile(item.path())) continue;
            try {
                const auto relative = relativePath(archiveRoot, item.path());
                if (relative.empty()) continue;
                present.insert(relative); ++result.discovered;
                const auto stamp = modifiedTicks(item.path());
                const auto bytes = item.file_size();
                const auto found = known.find(relative);
                if (found != known.end() && found->second.modifiedTicks == stamp &&
                    found->second.sizeBytes == bytes && found->second.summaryReady &&
                    found->second.summaryVersion == kArchiveSummaryVersion) {
                    ++result.unchanged; continue;
                }
                upsert(database, describe(archiveRoot, item.path())); ++result.indexed;
            } catch (...) {
                // A disappearing or corrupt archive is retried on the next reconcile.
                ++result.failed;
            }
        }
        sqlite3_stmt* remove{};
        sqlite3_prepare_v2(database.value,
            "DELETE FROM archive_objects WHERE relative_path=?", -1, &remove, nullptr);
        std::size_t removedArchiveObjects = 0;
        for (const auto& item : known) if (!present.count(item.first)) {
            sqlite3_reset(remove); sqlite3_clear_bindings(remove);
            sqlite3_bind_text(remove, 1, item.first.c_str(), -1, SQLITE_TRANSIENT);
            if (sqlite3_step(remove) == SQLITE_DONE) {
                ++result.removed; ++removedArchiveObjects;
            }
        }
        sqlite3_finalize(remove);
        result.removed += removeStaleSolutions(database, archiveRoot);
        if (result.indexed || removedArchiveObjects)
            invalidateCatalogSnapshots(database);
        database.execute("COMMIT");
    } catch (...) { database.execute("ROLLBACK"); throw; }
    return result;
}

}  // namespace dynsyslab::core
