#pragma once

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace dsl::compute {

enum class SolveNodeKind { batch, system, cluster, trajectory };

struct SolveNodeAddress {
    SolveNodeKind kind = SolveNodeKind::batch;
    std::size_t batch = 0;
    std::size_t system = 0;
    std::size_t cluster = 0;
    std::size_t trajectory = 0;
};

struct SolveSystemAddress {
    std::size_t batch = 0;
    std::size_t system = 0;
    friend bool operator==(const SolveSystemAddress& left,
                           const SolveSystemAddress& right) {
        return left.batch == right.batch && left.system == right.system;
    }
    friend bool operator!=(const SolveSystemAddress& left,
                           const SolveSystemAddress& right) {
        return !(left == right);
    }
};

struct SolveSystemDraft {
    std::string label{"New system"};
    std::string jobJson;
    std::size_t trajectories = 1;
    std::size_t clusters = 1;
};

struct SolveBatchDraft {
    std::string label;
    std::vector<SolveSystemDraft> systems;
};

[[nodiscard]] std::size_t assignedCluster(std::size_t trajectory,
                                          std::size_t trajectories,
                                          std::size_t clusters);

class SolveJobModel {
public:
    SolveJobModel();
    void reset(std::string systemLabel, std::string jobJson);
    std::size_t addBatch();
    SolveSystemAddress addSystem(std::size_t batch, std::string systemLabel = "New system");
    bool addCluster(const SolveSystemAddress& system);
    bool addTrajectory(const SolveSystemAddress& system);
    bool remove(const SolveNodeAddress& node);
    SolveSystemDraft* system(const SolveSystemAddress& address);
    const SolveSystemDraft* system(const SolveSystemAddress& address) const;
    const SolveBatchDraft* batch(std::size_t index) const;
    std::optional<SolveSystemAddress> owningSystem(const SolveNodeAddress& node) const;
    std::optional<SolveSystemAddress> firstSystem() const;
    const std::vector<SolveBatchDraft>& batches() const { return batches_; }
    std::vector<SolveBatchDraft>& batches() { return batches_; }

private:
    std::vector<SolveBatchDraft> batches_;
    std::string nextBatchLabel() const;
};

}  // namespace dsl::compute
