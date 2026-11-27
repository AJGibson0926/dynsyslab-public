#include <dynsyslab/core/CampaignRegistry.hpp>

#include "JobStoreInternal.hpp"

namespace dynsyslab::core {

using namespace jobstore;

struct CampaignRegistry::Implementation {
    sqlite3* database{};
    explicit Implementation(const std::filesystem::path& path) {
        std::filesystem::create_directories(path.parent_path());
        check(sqlite3_open_v2(path.u8string().c_str(), &database,
            SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE | SQLITE_OPEN_FULLMUTEX, nullptr),
            database, "open campaign registry");
        sqlite3_busy_timeout(database, 5000);
        execute(database, "PRAGMA journal_mode=WAL;");
        execute(database,
            "CREATE TABLE IF NOT EXISTS campaign_plans(id TEXT PRIMARY KEY,label TEXT NOT NULL,"
            "kind TEXT NOT NULL,point_count INTEGER NOT NULL,plan_json TEXT NOT NULL,"
            "created_utc TEXT NOT NULL);"
            "CREATE TABLE IF NOT EXISTS campaign_chunks(campaign_id TEXT NOT NULL,"
            "first_point INTEGER NOT NULL,point_count INTEGER NOT NULL,job_id TEXT NOT NULL DEFAULT '',"
            "manifest_path TEXT NOT NULL DEFAULT '',updated_utc TEXT NOT NULL,"
            "PRIMARY KEY(campaign_id,first_point));");
    }
    ~Implementation() { if (database) sqlite3_close(database); }
};

CampaignRegistry::CampaignRegistry(const std::filesystem::path& database)
    : implementation_(std::make_unique<Implementation>(database)) {}
CampaignRegistry::~CampaignRegistry() = default;

void CampaignRegistry::registerPlan(const CampaignPlanRecord& plan) {
    Statement(implementation_->database,
        "INSERT OR IGNORE INTO campaign_plans(id,label,kind,point_count,plan_json,created_utc)"
        " VALUES(?,?,?,?,?,?)", "register campaign plan")
        .text(plan.id).text(plan.label).text(plan.kind)
        .integer(static_cast<long long>(plan.pointCount)).text(plan.planJson)
        .text(plan.createdUtc.empty() ? timestamp() : plan.createdUtc).run();
}

void CampaignRegistry::recordChunk(const CampaignChunkRecord& chunk) {
    Statement(implementation_->database,
        "INSERT INTO campaign_chunks(campaign_id,first_point,point_count,job_id,manifest_path,"
        "updated_utc) VALUES(?,?,?,?,?,?) ON CONFLICT(campaign_id,first_point) DO UPDATE SET "
        "point_count=excluded.point_count,"
        "job_id=CASE WHEN excluded.job_id='' THEN job_id ELSE excluded.job_id END,"
        "manifest_path=CASE WHEN excluded.manifest_path='' THEN manifest_path"
        " ELSE excluded.manifest_path END,updated_utc=excluded.updated_utc", "record campaign chunk")
        .text(chunk.campaignId).integer(static_cast<long long>(chunk.firstPoint))
        .integer(static_cast<long long>(chunk.pointCount)).text(chunk.jobId)
        .text(chunk.manifestPath).text(timestamp()).run();
}

namespace {
CampaignPlanRecord planRecord(sqlite3_stmt* statement) {
    return {column(statement, 0), column(statement, 1), column(statement, 2),
            static_cast<std::size_t>(sqlite3_column_int64(statement, 3)), column(statement, 4),
            column(statement, 5)};
}
constexpr const char* planColumns = "id,label,kind,point_count,plan_json,created_utc";
}  // namespace

std::optional<CampaignPlanRecord> CampaignRegistry::plan(const std::string& id) const {
    Statement query(implementation_->database,
        std::string("SELECT ") + planColumns + " FROM campaign_plans WHERE id=?", "find campaign");
    query.text(id);
    if (!query.row()) return std::nullopt;
    return planRecord(query.get());
}

std::vector<CampaignPlanRecord> CampaignRegistry::plans(std::size_t limit) const {
    Statement query(implementation_->database,
        std::string("SELECT ") + planColumns
            + " FROM campaign_plans ORDER BY created_utc DESC,rowid DESC LIMIT ?", "list campaigns");
    query.integer(static_cast<long long>(limit));
    std::vector<CampaignPlanRecord> result;
    while (query.row()) result.push_back(planRecord(query.get()));
    return result;
}

std::vector<CampaignChunkRecord> CampaignRegistry::chunks(const std::string& id) const {
    Statement query(implementation_->database,
        "SELECT campaign_id,first_point,point_count,job_id,manifest_path,updated_utc"
        " FROM campaign_chunks WHERE campaign_id=? ORDER BY first_point", "list campaign chunks");
    query.text(id);
    std::vector<CampaignChunkRecord> result;
    while (query.row())
        result.push_back({column(query.get(), 0),
            static_cast<std::size_t>(sqlite3_column_int64(query.get(), 1)),
            static_cast<std::size_t>(sqlite3_column_int64(query.get(), 2)),
            column(query.get(), 3), column(query.get(), 4), column(query.get(), 5)});
    return result;
}

}  // namespace dynsyslab::core
