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

// Read-only, bounded archive inspection. It never rewrites or removes data.
[[nodiscard]] ArchiveAuditReport auditArchives(
    const std::filesystem::path& root,
    const ArchiveAuditScope& scope = {});

// Stable machine-readable representation used by the compiled CLI and UI.
[[nodiscard]] std::string archiveAuditJson(const ArchiveAuditReport& report,
                                           const ArchiveAuditScope& scope = {});

}  // namespace dynsyslab::core
