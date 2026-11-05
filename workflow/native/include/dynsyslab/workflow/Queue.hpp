#pragma once

#include <dynsyslab/core/Execution.hpp>
#include <dynsyslab/core/JobStore.hpp>
#include <dynsyslab/workflow/Json.hpp>
#include <dynsyslab/workflow/Run.hpp>

#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace dynsyslab::workflow {

struct QueueRunResult {
    enum class State { empty, completed, cancelled };
    State state{State::empty};
    std::optional<dynsyslab::core::JobRecord> job;
    std::optional<RunResult> run;
};

class JobQueue {
public:
    explicit JobQueue(std::filesystem::path archiveRoot);

    [[nodiscard]] std::string submit(const json::Value& job);
    [[nodiscard]] std::vector<dynsyslab::core::JobRecord> list(
        std::size_t limit = 100) const;
    [[nodiscard]] std::optional<dynsyslab::core::JobRecord> find(
        const std::string& id) const;
    [[nodiscard]] dynsyslab::core::JobRecord cancel(const std::string& id);
    [[nodiscard]] std::size_t recoverStale(double ageSeconds);
    [[nodiscard]] QueueRunResult runNext(
        const std::string& worker,
        const std::filesystem::path& progressPath = {},
        dynsyslab::core::ExecutionLimits limits = {},
        dynsyslab::core::CancellationProbe externalCancellation = {});

private:
    std::filesystem::path archiveRoot_;
    std::filesystem::path database_;
};

[[nodiscard]] json::Value queueRunJson(const RunResult& run);

}  // namespace dynsyslab::workflow
