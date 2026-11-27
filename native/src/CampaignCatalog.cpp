// Campaign catalog rebuild (core/CampaignCatalog.hpp).
#include <dynsyslab/core/CampaignCatalog.hpp>

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/ArchiveAudit.hpp>
#include <dynsyslab/core/Json.hpp>

#include <map>
#include <set>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

struct Campaign {
    json::Value::Object summary;
    std::set<std::string> manifests, archives, missing;
};

json::Value strings(const std::set<std::string>& values) {
    json::Value::Array result;
    for (const auto& value : values) result.emplace_back(value);
    return json::Value(std::move(result));
}

bool isFile(const std::filesystem::path& path) {
    std::error_code ignored;
    return std::filesystem::is_regular_file(path, ignored);
}

}  // namespace

std::string campaignCatalogJson(const std::filesystem::path& archiveRoot) {
    const auto root = std::filesystem::absolute(archiveRoot).lexically_normal();
    std::map<std::string, Campaign> campaigns;
    json::Value::Array unreadable;
    std::error_code ignored;
    if (std::filesystem::is_directory(root / L"campaigns", ignored))
        for (const auto& entry : std::filesystem::directory_iterator(root / L"campaigns")) {
            if (!entry.is_regular_file() || entry.path().extension() != ".json") continue;
            try {
                const auto manifest = json::load(entry.path().wstring());
                const auto id = manifest.at("id").string(manifest.at("campaign_id").string());
                if (id.empty()) continue;
                auto& campaign = campaigns[id];
                for (const char* key : {"label", "kind", "created_utc", "schema"})
                    if (!manifest.at(key).isNull()) campaign.summary[key] = manifest.at(key);
                campaign.summary["points"] = json::Value(static_cast<double>(manifest.at("points").array().size()));
                campaign.manifests.insert(entry.path().generic_u8string());
                for (const auto& point : manifest.at("points").array())
                    for (const auto& archive : point.at("solution_archives").array())
                        if (!isFile(std::filesystem::u8path(archive.string())))
                            campaign.missing.insert(archive.string());
            } catch (const std::exception& error) {
                unreadable.emplace_back(json::Value::Object{{"path", json::Value(entry.path().generic_u8string())},
                                                            {"error", json::Value(error.what())}});
            }
        }
    ArchiveAuditScope scope;
    scope.modifiedWithinDays = -1.0;
    scope.sections = {"solutions"};
    scope.includeDerived = true;
    for (const auto& path : scopedArchives(root, scope)) {
        try {
            const auto text = readNpzIndexStored(path).metadataJson;
            const auto id = text.empty() ? std::string{} : json::parse(text).at("campaign_id").string();
            if (!id.empty()) campaigns[id].archives.insert(path.generic_u8string());
        } catch (const std::exception& error) {
            unreadable.emplace_back(json::Value::Object{{"path", json::Value(path.generic_u8string())},
                                                        {"error", json::Value(error.what())}});
        }
    }
    json::Value::Array list;
    std::size_t archives = 0;
    for (auto& [id, campaign] : campaigns) {
        auto summary = std::move(campaign.summary);
        summary["id"] = json::Value(id);
        summary["manifests"] = strings(campaign.manifests);
        summary["solution_archives"] = strings(campaign.archives);
        summary["missing_archives"] = strings(campaign.missing);
        archives += campaign.archives.size();
        list.emplace_back(std::move(summary));
    }
    return json::dump(json::Value(json::Value::Object{{"schema", json::Value(kCampaignCatalogSchema)},
        {"root", json::Value(root.generic_u8string())},
        {"campaigns", json::Value(std::move(list))},
        {"campaign_count", json::Value(static_cast<double>(campaigns.size()))},
        {"solution_archive_count", json::Value(static_cast<double>(archives))},
        {"unreadable", json::Value(std::move(unreadable))}}));
}

}  // namespace dynsyslab::core
