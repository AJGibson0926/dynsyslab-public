#pragma once

#include "dynsyslab_compute/Registry.hpp"
#include "dynsyslab_viewer/Model.hpp"
#include <dynsyslab/core/ArchiveAudit.hpp>
#include <dynsyslab/core/LegacyCsv.hpp>
#include <dynsyslab/workflow/Queue.hpp>

#include <filesystem>
#include <functional>
#include <string>
#include <vector>

namespace dsl::compute {

struct BackendConfig {
    std::filesystem::path projectRoot;
    std::filesystem::path archiveRoot;
    std::filesystem::path registry;
    std::filesystem::path catalog;
    std::filesystem::path cache;
};

struct ProgressStatus {
    bool available = false;
    bool determinate = false;
    double percent = 0.0;
    double estimatedRemainingSeconds = -1.0;
    int completedTrajectories = 0;
    int plannedTrajectories = 0;
    std::wstring state;
    std::wstring message;
};

class Backend {
public:
    explicit Backend(BackendConfig config) : config_(std::move(config)) {}
    Registry loadRegistry() const;
    Registry loadRegistry(const dsl::Catalog& catalog) const;
    dsl::Catalog loadCatalog() const;
    dsl::Catalog refreshCatalog() const;
    std::wstring execute(const std::string& action, const std::string& jobJson,
                         const std::filesystem::path& progress = {},
                         dynsyslab::core::CancellationProbe cancellation = {}) const;
    std::vector<dynsyslab::core::JobRecord> queueJobs() const;
    std::string submitQueued(const std::string& jobJson) const;
    dynsyslab::core::JobRecord cancelQueued(const std::string& id) const;
    std::size_t recoverQueued(double ageSeconds) const;
    dynsyslab::workflow::QueueRunResult runNextQueued(
        const std::filesystem::path& progress = {},
        dynsyslab::core::CancellationProbe cancellation = {}) const;
    std::wstring auditArchive(const dynsyslab::core::ArchiveAuditScope& scope) const;
    std::wstring deduplicateArchive(
        const dynsyslab::core::ArchiveAuditScope& scope, bool apply) const;
    std::wstring importLegacyCsv(
        const std::filesystem::path& source,
        const dynsyslab::core::LegacyCsvImportOptions& options) const;
    ProgressStatus progress(const std::filesystem::path& path) const;
    bool openViewer() const;
    const BackendConfig& config() const { return config_; }
private:
    BackendConfig config_;
};

}  // namespace dsl::compute
