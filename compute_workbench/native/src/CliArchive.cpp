#include "CliArchive.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/ArchiveAudit.hpp>
#include <dynsyslab/core/ArchiveDeduplication.hpp>
#include <dynsyslab/core/ArchiveIndex.hpp>
#include <dynsyslab/core/ArchiveMaintenance.hpp>
#include <dynsyslab/core/CampaignCatalog.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/OutputPath.hpp>
#include <dynsyslab/core/StorageLocations.hpp>
#include <dynsyslab/core/Text.hpp>

#include <fstream>
#include <iostream>
#include <stdexcept>

namespace dynsyslab::cli {
namespace {
namespace json = dynsyslab::core::json;

core::ArchiveAuditScope scopeOf(const Arguments& args) {
    core::ArchiveAuditScope scope;
    scope.modifiedWithinDays = args.auditDays;
    if (!args.auditSections.empty()) scope.sections = args.auditSections;
    scope.systems = args.auditSystems;
    scope.campaignId = args.auditCampaign;
    scope.includeDerived = args.auditDerived;
    return scope;
}

std::filesystem::path positional(const Arguments& args, std::size_t index, const char* usage) {
    if (args.positional.size() <= index) throw std::invalid_argument(usage);
    return std::filesystem::path(args.positional[index]);
}

void indexUpdate(const Arguments& args) {
    const auto update = core::reconcileArchiveIndex(args.archiveRoot);
    json::Value::Object result;
    result["schema"] = "dynsyslab.archive-index-update/v1";
    result["index"] = core::archiveIndexPath(args.archiveRoot).generic_u8string();
    result["discovered"] = static_cast<double>(update.discovered);
    result["unchanged"] = static_cast<double>(update.unchanged);
    result["indexed"] = static_cast<double>(update.indexed);
    result["removed"] = static_cast<double>(update.removed);
    result["failed"] = static_cast<double>(update.failed);
    std::cout << json::dump(json::Value(std::move(result))) << '\n';
}

void audit(const Arguments& args) {
    const auto scope = scopeOf(args);
    const auto report = core::auditArchives(args.archiveRoot, scope);
    if (args.report.empty()) {
        std::cout << core::archiveAuditJson(report, scope) << '\n';
        return;
    }
    const auto files = core::writeArchiveAuditReport(report, scope, args.report);
    auto result = json::parse(core::archiveAuditJson(report, scope));
    result.at("report_json") = json::Value(files.json.generic_u8string());
    result.at("report_markdown") = json::Value(files.markdown.generic_u8string());
    std::cout << json::dump(result) << '\n';
}

std::vector<std::filesystem::path> selectedFiles(const Arguments& args) {
    if (args.positional.size() <= 2) return core::scopedArchives(args.archiveRoot, scopeOf(args));
    std::vector<std::filesystem::path> files;
    for (std::size_t index = 2; index < args.positional.size(); ++index)
        files.emplace_back(core::resolveOutputPath(std::filesystem::path(args.positional[index])));
    return files;
}

void editMetadata(const Arguments& args) {
    constexpr const char* usage = "archive edit-metadata needs <source.npz> <output.npz> --metadata-json JSON";
    const auto source = positional(args, 2, usage);
    const auto output = positional(args, 3, usage);
    if (args.metadataPatch.empty()) throw std::invalid_argument(usage);
    const auto written = core::editArchiveMetadata(source, output, args.metadataPatch);
    (void)core::indexArchiveFiles(args.archiveRoot, {written});
    std::cout << json::dump(json::Value(json::Value::Object{{"schema", json::Value("dynsyslab.archive-metadata-edit/v1")},
        {"source", json::Value(std::filesystem::absolute(source).generic_u8string())},
        {"output", json::Value(written.generic_u8string())}})) << '\n';
}

void campaigns(const Arguments& args) {
    const auto text = core::campaignCatalogJson(args.archiveRoot);
    if (!args.output.empty()) {
        const auto path = core::prepareOutputFile(args.output, "--output");
        std::ofstream stream(path, std::ios::binary | std::ios::trunc);
        if (!stream || !(stream << text << '\n')) throw std::runtime_error("unable to write " + path.u8string());
    }
    std::cout << text << '\n';
}

}  // namespace

void archiveCommand(const Arguments& args) {
    constexpr const char* usage =
        "archive requires index, audit, deduplicate, recompress, share-time-grids, edit-metadata or campaigns";
    if (args.positional.size() < 2) throw std::invalid_argument(usage);
    const auto& operation = args.positional[1];
    if (operation == L"index") indexUpdate(args);
    else if (operation == L"audit") audit(args);
    else if (operation == L"deduplicate") {
        const auto result = core::deduplicateArchives(args.archiveRoot, scopeOf(args), args.apply);
        if (args.apply) (void)core::reconcileArchiveIndex(args.archiveRoot);
        std::cout << core::archiveDeduplicationJson(result) << '\n';
    } else if (operation == L"recompress") {
        // Lossless re-encoding; LZMA unless --compression names another codec.
        const auto codec = args.compression.empty() ? core::ArchiveCompression::lzma
                                                    : core::parseArchiveCompression(args.compression);
        const auto result = core::recompressArchives(selectedFiles(args), codec, args.apply);
        if (args.apply) (void)core::reconcileArchiveIndex(args.archiveRoot);
        std::cout << core::archiveMaintenanceJson(result) << '\n';
    } else if (operation == L"share-time-grids") {
        const auto result = core::shareTimeGrids(args.archiveRoot, selectedFiles(args), args.apply);
        if (args.apply) (void)core::reconcileArchiveIndex(args.archiveRoot);
        std::cout << core::archiveMaintenanceJson(result) << '\n';
    } else if (operation == L"edit-metadata") editMetadata(args);
    else if (operation == L"campaigns") campaigns(args);
    else throw std::invalid_argument(usage);
}

void storageCommand(const Arguments& args) {
    constexpr const char* usage = "storage requires list, register, copy, relocate or append";
    if (args.positional.size() < 2) throw std::invalid_argument(usage);
    const auto& operation = args.positional[1];
    const auto& root = args.archiveRoot;
    if (operation == L"list") {
        std::cout << core::storageLocationsJson(core::storageLocations(root)) << '\n';
    } else if (operation == L"register") {
        constexpr const char* need = "storage register needs <name> <folder>";
        core::StorageLocation location{core::narrow(positional(args, 2, need).wstring()),
            core::resolveOutputPath(positional(args, 3, need)), !args.readOnly, args.purpose};
        core::registerStorageLocation(root, location);
        std::cout << core::storageLocationsJson(core::storageLocations(root)) << '\n';
    } else if (operation == L"copy") {
        constexpr const char* need = "storage copy needs <file> <destination folder>";
        const auto copy = core::copyVerified(core::resolveOutputPath(positional(args, 2, need)),
                                             core::resolveOutputPath(positional(args, 3, need)), args.removeSource);
        std::cout << core::verifiedCopyJson(copy) << '\n';
    } else if (operation == L"relocate") {
        constexpr const char* need = "storage relocate needs <campaign-id> <registered location>";
        const auto campaign = core::narrow(positional(args, 2, need).wstring());
        const auto destination = core::narrow(positional(args, 3, need).wstring());
        const auto files = core::campaignArchives(root, campaign);
        const auto plan = core::createRelocationPlan(root, campaign, destination);
        json::Value::Array copies;
        for (const auto& file : files) copies.push_back(json::parse(core::verifiedCopyJson(core::appendRelocationCopy(plan, file))));
        std::cout << json::dump(json::Value(json::Value::Object{{"schema", json::Value("dynsyslab.campaign-relocation-result/v1")},
            {"plan", json::Value(plan.generic_u8string())}, {"campaign_id", json::Value(campaign)},
            {"copies", json::Value(std::move(copies))}})) << '\n';
    } else if (operation == L"append") {
        constexpr const char* need = "storage append needs <plan.json> <file>...";
        const auto plan = core::resolveOutputPath(positional(args, 2, need));
        (void)positional(args, 3, need);
        json::Value::Array copies;
        for (std::size_t index = 3; index < args.positional.size(); ++index)
            copies.push_back(json::parse(core::verifiedCopyJson(core::appendRelocationCopy(
                plan, core::resolveOutputPath(std::filesystem::path(args.positional[index]))))));
        std::cout << json::dump(json::Value(json::Value::Object{{"schema", json::Value("dynsyslab.campaign-relocation-result/v1")},
            {"plan", json::Value(plan.generic_u8string())}, {"copies", json::Value(std::move(copies))}})) << '\n';
    } else throw std::invalid_argument(usage);
}

}  // namespace dynsyslab::cli
