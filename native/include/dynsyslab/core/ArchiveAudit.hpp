#pragma once

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace dynsyslab::core {

struct ArchiveAuditScope {
    // Negative means no age limit. The default mirrors the conservative
    // interactive audit window used by the retained client.
    double modifiedWithinDays{2.0};
    std::vector<std::string> sections{"solutions"};
    std::vector<std::string> systems;
    // Non-empty: only archives whose metadata campaign_id equals it (written by native campaigns).
    std::string campaignId;
    bool includeDerived{};
};

struct ArchiveAuditIssue {
    std::filesystem::path path;
    std::string message;
};

struct ArchiveDuplicateGroup {
    std::string kind;
    std::vector<std::filesystem::path> files;
    std::uintmax_t redundantBytes{};
};

struct ArchiveAuditReport {
    std::filesystem::path root;
    std::size_t scannedFiles{};
    std::uintmax_t totalBytes{};
    std::vector<ArchiveDuplicateGroup> exactDuplicateGroups;
    std::vector<ArchiveDuplicateGroup> samePlanGroups;
    std::vector<ArchiveAuditIssue> unreadableFiles;
};

// The completed .npz files the scope selects, sorted; files whose metadata cannot be read for
// the campaign filter go to `unreadable` when given. Shared by audit, deduplication and the
// storage maintenance commands (ArchiveMaintenance.hpp).
[[nodiscard]] std::vector<std::filesystem::path> scopedArchives(
    const std::filesystem::path& root, const ArchiveAuditScope& scope,
    std::vector<ArchiveAuditIssue>* unreadable = nullptr);

// Read-only, bounded archive inspection. It never rewrites or removes data.
[[nodiscard]] ArchiveAuditReport auditArchives(
    const std::filesystem::path& root,
    const ArchiveAuditScope& scope = {});

// Stable machine-readable representation used by the compiled CLI and UI.
[[nodiscard]] std::string archiveAuditJson(const ArchiveAuditReport& report,
                                           const ArchiveAuditScope& scope = {});
// Human-readable summary of the same report.
[[nodiscard]] std::string archiveAuditMarkdown(const ArchiveAuditReport& report,
                                               const ArchiveAuditScope& scope = {});

struct ArchiveAuditReportFiles {
    std::filesystem::path json;
    std::filesystem::path markdown;
};
// Writes archive-audit__<UTC stamp>.json and .md into `directory` (created when missing).
ArchiveAuditReportFiles writeArchiveAuditReport(const ArchiveAuditReport& report,
                                                const ArchiveAuditScope& scope,
                                                const std::filesystem::path& directory);

}  // namespace dynsyslab::core
