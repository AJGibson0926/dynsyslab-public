#include "dynsyslab_compute/SolveJobModel.hpp"

#include <algorithm>

namespace dsl::compute {

std::size_t assignedCluster(std::size_t trajectory, std::size_t trajectories,
                            std::size_t clusters) {
    if (!trajectories || !clusters) return 0;
    return std::min(clusters - 1, trajectory * clusters / trajectories);
}

SolveJobModel::SolveJobModel() { batches_.push_back({"Batch 1", {}}); }

void SolveJobModel::reset(std::string systemLabel, std::string jobJson) {
    batches_.clear();
    SolveBatchDraft batch{"Batch 1", {}};
    batch.systems.push_back({std::move(systemLabel), std::move(jobJson), 1, 1});
    batches_.push_back(std::move(batch));
}

std::string SolveJobModel::nextBatchLabel() const {
    return "Batch " + std::to_string(batches_.size() + 1);
}

std::size_t SolveJobModel::addBatch() {
    batches_.push_back({nextBatchLabel(), {}});
    return batches_.size() - 1;
}

SolveSystemAddress SolveJobModel::addSystem(std::size_t batch, std::string systemLabel) {
    if (batch >= batches_.size()) batch = batches_.empty() ? addBatch() : batches_.size() - 1;
    auto& systems = batches_[batch].systems;
    systems.push_back({std::move(systemLabel), {}, 1, 1});
    return {batch, systems.size() - 1};
}

SolveSystemDraft* SolveJobModel::system(const SolveSystemAddress& address) {
    if (address.batch >= batches_.size()
        || address.system >= batches_[address.batch].systems.size()) return nullptr;
    return &batches_[address.batch].systems[address.system];
}

const SolveSystemDraft* SolveJobModel::system(const SolveSystemAddress& address) const {
    if (address.batch >= batches_.size()
        || address.system >= batches_[address.batch].systems.size()) return nullptr;
    return &batches_[address.batch].systems[address.system];
}

const SolveBatchDraft* SolveJobModel::batch(std::size_t index) const {
    return index < batches_.size() ? &batches_[index] : nullptr;
}

bool SolveJobModel::addCluster(const SolveSystemAddress& address) {
    auto* draft = system(address);
    if (!draft) return false;
    ++draft->clusters;
    draft->trajectories = std::max(draft->trajectories, draft->clusters);
    return true;
}

bool SolveJobModel::addTrajectory(const SolveSystemAddress& address) {
    auto* draft = system(address);
    if (!draft) return false;
    ++draft->trajectories;
    return true;
}

std::optional<SolveSystemAddress> SolveJobModel::owningSystem(
    const SolveNodeAddress& node) const {
    if (node.kind == SolveNodeKind::batch) return std::nullopt;
    const SolveSystemAddress result{node.batch, node.system};
    return system(result) ? std::optional<SolveSystemAddress>(result) : std::nullopt;
}

std::optional<SolveSystemAddress> SolveJobModel::firstSystem() const {
    for (std::size_t batch = 0; batch < batches_.size(); ++batch)
        if (!batches_[batch].systems.empty()) return SolveSystemAddress{batch, 0};
    return std::nullopt;
}

bool SolveJobModel::remove(const SolveNodeAddress& node) {
    if (node.batch >= batches_.size()) return false;
    if (node.kind == SolveNodeKind::batch) {
        if (batches_.size() == 1) return false;
        batches_.erase(batches_.begin() + static_cast<std::ptrdiff_t>(node.batch));
        return true;
    }
    const SolveSystemAddress address{node.batch, node.system};
    auto* draft = system(address);
    if (!draft) return false;
    if (node.kind == SolveNodeKind::system) {
        auto& systems = batches_[node.batch].systems;
        systems.erase(systems.begin() + static_cast<std::ptrdiff_t>(node.system));
        return true;
    }
    if (node.kind == SolveNodeKind::cluster) {
        if (draft->clusters == 1) return false;
        --draft->clusters;
        return true;
    }
    if (draft->trajectories == 1) return false;
    --draft->trajectories;
    draft->clusters = std::min(draft->clusters, draft->trajectories);
    return true;
}

}  // namespace dsl::compute
