#pragma once

#include <cstddef>
#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace dynsyslab::core {

enum class JobStatus { queued, running, completed, failed, cancelled };

struct JobRecord {
    std::string id;
    std::string label;
    std::string kind;
    std::string manifestJson;
    JobStatus status{JobStatus::queued};
    double progress{};
    std::string message;
    std::string worker;
    std::string createdUtc;
    std::string updatedUtc;
    std::size_t attempts{};
};

class JobStore {
public:
    explicit JobStore(const std::filesystem::path& path);
    ~JobStore();
    JobStore(const JobStore&) = delete;
    JobStore& operator=(const JobStore&) = delete;

    [[nodiscard]] std::string submit(const std::string& label, const std::string& kind,
                                     const std::string& manifestJson);
    [[nodiscard]] std::optional<JobRecord> claimNext(const std::string& worker);
    void updateProgress(const std::string& id, double progress, const std::string& message);
    void complete(const std::string& id, const std::string& message);
    void fail(const std::string& id, const std::string& message);
    void cancel(const std::string& id, const std::string& message = {});
    [[nodiscard]] std::size_t recoverStale(double ageSeconds);
    [[nodiscard]] std::vector<JobRecord> list(
        std::size_t limit = 100, std::optional<JobStatus> status = std::nullopt) const;
    [[nodiscard]] std::optional<JobRecord> find(const std::string& id) const;

private:
    struct Implementation;
    std::unique_ptr<Implementation> implementation_;
};

[[nodiscard]] std::string jobStatusName(JobStatus status);

}  // namespace dynsyslab::core
