#include "dynsyslab_compute/Backend.hpp"
#include <dynsyslab/core/ArchiveDeduplication.hpp>
#include <dynsyslab/core/ArchiveIndex.hpp>
#include <dynsyslab/workflow/Plan.hpp>
#include <dynsyslab/workflow/Run.hpp>
#include "dynsyslab_viewer/Json.hpp"
#include "dynsyslab_viewer/Process.hpp"

#include <windows.h>

#include <stdexcept>
#include <utility>

namespace dsl::compute {

namespace {
constexpr wchar_t kViewerWindowClass[] = L"DynSysLabViewerWorkbench";

std::filesystem::path viewerExecutable(const BackendConfig& config) {
    wchar_t buffer[32768]{};
    const DWORD count = GetModuleFileNameW(nullptr, buffer, 32768);
    const auto sibling = std::filesystem::path(std::wstring(buffer, count)).parent_path()
        / L"DynSysLabViewer.exe";
    if (std::filesystem::is_regular_file(sibling)) return sibling;
    const auto development = config.projectRoot / L"build" / L"viewer-native"
        / L"Release" / L"DynSysLabViewer.exe";
    return std::filesystem::is_regular_file(development) ? development : sibling;
}

bool activateExistingViewer() {
    const HWND viewer = FindWindowW(kViewerWindowClass, nullptr);
    if (!viewer) return false;
    if (IsIconic(viewer)) ShowWindow(viewer, SW_RESTORE);
    else ShowWindow(viewer, SW_SHOW);
    BringWindowToTop(viewer);
    SetForegroundWindow(viewer);
    return true;
}
}

Registry Backend::loadRegistry() const {
    return loadRegistry(loadCatalog());
}

Registry Backend::loadRegistry(const dsl::Catalog& catalog) const {
    if (!config_.registry.empty() && std::filesystem::is_regular_file(config_.registry)) {
        try { return Registry::load(config_.registry.wstring()); }
        catch (const std::exception&) { /* Prefer a valid built-in native registry. */ }
    }
    return Registry::builtin(catalog);
}
dsl::Catalog Backend::loadCatalog() const {
    auto indexed = dsl::Catalog::indexed(config_.projectRoot, config_.archiveRoot);
    if (!indexed.items.empty()) return indexed;
    if (!config_.catalog.empty() && std::filesystem::is_regular_file(config_.catalog)) {
        try { return dsl::Catalog::load(config_.catalog.wstring()); }
        catch (const std::exception&) { /* Fall through to native discovery. */ }
    }
    return dsl::Catalog::scan(config_.projectRoot, config_.archiveRoot);
}

dsl::Catalog Backend::refreshCatalog() const {
    (void)dynsyslab::core::reconcileArchiveIndex(config_.archiveRoot);
    auto indexed = dsl::Catalog::indexed(config_.projectRoot, config_.archiveRoot);
    return indexed.items.empty()
        ? dsl::Catalog::scan(config_.projectRoot, config_.archiveRoot)
        : indexed;
}

std::wstring Backend::execute(const std::string& action, const std::string& jobJson,
                              const std::filesystem::path& progress,
                              dynsyslab::core::CancellationProbe cancellation) const {
    const auto root = dsl::json::parse(jobJson);
    dynsyslab::workflow::validateJob(root);
    const auto kind = root.at("kind").string();
    if (action == "plan") return L"Native plan\r\n\r\n"
        + dsl::wide(dsl::json::dump(dynsyslab::workflow::planJob(root)))
        + L"\r\n\r\nPython required: no";
    dynsyslab::core::ExecutionContext execution;execution.cancellation=std::move(cancellation);
    return dynsyslab::workflow::executeJob(root,config_.archiveRoot,progress,
        execution.cancellation?&execution:nullptr).summary;
}

std::vector<dynsyslab::core::JobRecord> Backend::queueJobs() const {
    return dynsyslab::workflow::JobQueue(config_.archiveRoot).list();
}

std::string Backend::submitQueued(const std::string& jobJson) const {
    return dynsyslab::workflow::JobQueue(config_.archiveRoot).submit(
        dsl::json::parse(jobJson));
}

dynsyslab::core::JobRecord Backend::cancelQueued(const std::string& id) const {
    return dynsyslab::workflow::JobQueue(config_.archiveRoot).cancel(id);
}

std::size_t Backend::recoverQueued(double ageSeconds) const {
    return dynsyslab::workflow::JobQueue(config_.archiveRoot).recoverStale(ageSeconds);
}

dynsyslab::workflow::QueueRunResult Backend::runNextQueued(
    const std::filesystem::path& progress,
    dynsyslab::core::CancellationProbe cancellation) const {
    return dynsyslab::workflow::JobQueue(config_.archiveRoot).runNext(
        "dynsyslab-native-ui",progress,{},std::move(cancellation));
}

std::wstring Backend::auditArchive(
    const dynsyslab::core::ArchiveAuditScope& scope) const {
    return dsl::wide(dynsyslab::core::archiveAuditJson(
        dynsyslab::core::auditArchives(config_.archiveRoot, scope), scope));
}

std::wstring Backend::deduplicateArchive(
    const dynsyslab::core::ArchiveAuditScope& scope, bool apply) const {
    const auto result = dynsyslab::core::deduplicateArchives(
        config_.archiveRoot, scope, apply);
    if (apply) (void)dynsyslab::core::reconcileArchiveIndex(config_.archiveRoot);
    return dsl::wide(dynsyslab::core::archiveDeduplicationJson(result));
}

std::wstring Backend::importLegacyCsv(
    const std::filesystem::path& source,
    const dynsyslab::core::LegacyCsvImportOptions& options) const {
    const auto result = dynsyslab::core::importLegacyDynBatchCsv(
        source, config_.archiveRoot, options);
    (void)dynsyslab::core::indexArchiveFiles(config_.archiveRoot, result.archives);
    return dsl::wide(dsl::json::dump(result.summary()));
}

ProgressStatus Backend::progress(const std::filesystem::path& path) const {
    ProgressStatus result;
    if (path.empty() || !std::filesystem::is_regular_file(path)) return result;
    try {
        const auto root = dsl::json::load(path.wstring());
        if (root.at("schema").string() != "dynsyslab.compute-progress/v1") return result;
        result.available = true;
        result.determinate = root.at("determinate").boolean();
        result.percent = root.at("percent").number();
        result.estimatedRemainingSeconds = root.at("estimated_remaining_seconds").number(-1.0);
        result.completedTrajectories = static_cast<int>(root.at("completed_trajectories").number());
        result.plannedTrajectories = static_cast<int>(root.at("planned_trajectories").number());
        result.state = dsl::wide(root.at("state").string());
        result.message = dsl::wide(root.at("message").string());
    } catch (...) {
        // Atomic writers make this rare; a later UI poll retries naturally.
    }
    return result;
}

bool Backend::openViewer() const {
    if (activateExistingViewer()) return true;
    const auto executable = viewerExecutable(config_);
    const auto viewerCache = config_.cache.parent_path() / L"Viewer";
    std::vector<std::wstring> command{executable.wstring(), L"--project-root",
        config_.projectRoot.wstring(), L"--archive-root", config_.archiveRoot.wstring(),
        L"--cache", viewerCache.wstring()};
    return dsl::launchDetached(command, config_.projectRoot.wstring()).exitCode == 0;
}

}  // namespace dsl::compute
