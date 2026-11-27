// Audit report files: archive-audit__<stamp>.json (archiveAuditJson) and a Markdown summary.
#include <dynsyslab/core/ArchiveAudit.hpp>

#include <dynsyslab/core/OutputPath.hpp>

#include <chrono>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

std::string utcStamp() {
    const std::time_t raw = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::tm utc{};
#ifdef _WIN32
    gmtime_s(&utc, &raw);
#else
    gmtime_r(&raw, &utc);
#endif
    std::ostringstream out;
    out << std::put_time(&utc, "%Y%m%dT%H%M%SZ");
    return out.str();
}

void writeText(const std::filesystem::path& path, const std::string& text) {
    const auto temporary = path.parent_path() / (path.filename().wstring() + L".partial");
    {
        std::ofstream stream(temporary, std::ios::binary | std::ios::trunc);
        if (!stream || !(stream << text) || !stream.flush())
            throw std::runtime_error("unable to write audit report " + path.u8string());
    }
    std::filesystem::rename(temporary, path);
}

std::string join(const std::vector<std::string>& values) {
    std::string text;
    for (const auto& value : values) text += (text.empty() ? "" : ", ") + value;
    return text.empty() ? "all" : text;
}

}  // namespace

std::string archiveAuditMarkdown(const ArchiveAuditReport& report, const ArchiveAuditScope& scope) {
    std::uintmax_t redundantBytes{}; std::size_t redundantFiles{};
    for (const auto& group : report.exactDuplicateGroups) {
        redundantBytes += group.redundantBytes; redundantFiles += group.files.size() - 1;
    }
    std::ostringstream out;
    out << "# DynSysLab archive audit\n\n"
        << "- Root: `" << report.root.generic_u8string() << "`\n"
        << "- Scope: sections " << join(scope.sections) << "; systems " << join(scope.systems)
        << "; campaign " << (scope.campaignId.empty() ? "all" : scope.campaignId) << "; age "
        << (scope.modifiedWithinDays < 0.0 ? std::string("all")
                                            : "within " + std::to_string(scope.modifiedWithinDays) + " days")
        << (scope.includeDerived ? "; derived included" : "; derived excluded") << "\n"
        << "- Scanned files: " << report.scannedFiles << " (" << report.totalBytes << " bytes)\n"
        << "- Exact duplicate groups: " << report.exactDuplicateGroups.size() << "\n"
        << "- Redundant files: " << redundantFiles << "\n"
        << "- Reclaimable storage: " << std::fixed << std::setprecision(2)
        << static_cast<double>(redundantBytes) / 1048576.0 << " MiB\n\n";
    if (!report.exactDuplicateGroups.empty()) {
        out << "## Exact duplicate solutions\n\nIdentical `time` and `values` arrays.\n\n";
        std::size_t index = 0;
        for (const auto& group : report.exactDuplicateGroups) {
            out << "### Group " << ++index << ": " << group.files.size() << " files\n\n";
            for (const auto& file : group.files) out << "- `" << file.generic_u8string() << "`\n";
            out << '\n';
        }
    }
    if (!report.samePlanGroups.empty())
        out << "## Repeated plans\n\n" << report.samePlanGroups.size()
            << " groups share a numerical plan (system, parameters, initial state, solver, forcing);"
               " the exact groups above are a subset.\n\n";
    if (!report.unreadableFiles.empty()) {
        out << "## Skipped files\n\n";
        for (const auto& issue : report.unreadableFiles)
            out << "- `" << issue.path.generic_u8string() << "`: " << issue.message << '\n';
        out << '\n';
    }
    out << "The audit is read-only. `dynsyslab-cli archive deduplicate --apply` moves unreferenced exact"
           " duplicates to archive trash.\n";
    return out.str();
}

ArchiveAuditReportFiles writeArchiveAuditReport(const ArchiveAuditReport& report,
                                                const ArchiveAuditScope& scope,
                                                const std::filesystem::path& directory) {
    const auto folder = prepareOutputDirectory(directory, "--report");
    const auto base = "archive-audit__" + utcStamp();
    ArchiveAuditReportFiles files{folder / std::filesystem::u8path(base + ".json"),
                                  folder / std::filesystem::u8path(base + ".md")};
    writeText(files.json, archiveAuditJson(report, scope) + "\n");
    writeText(files.markdown, archiveAuditMarkdown(report, scope));
    return files;
}

}  // namespace dynsyslab::core
