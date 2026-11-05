#pragma once

#include "dynsyslab_compute/SolveJobModel.hpp"
#include "dynsyslab_viewer/Theme.hpp"

#include <windows.h>

#include <optional>
#include <string>
#include <vector>

namespace dsl::compute {

constexpr int kSolveTreeId = 5100;
constexpr int kAddBatchId = 5101;
constexpr int kAddSystemId = 5102;
constexpr int kAddClusterId = 5103;
constexpr int kAddTrajectoryId = 5104;
constexpr int kRemoveSolveNodeId = 5105;
constexpr int kSolveJobNameId = 5106;

class SolveJobPane {
public:
    bool create(HWND parent, const dsl::theme::Fonts& fonts,
                int campaignCommandId);
    void show(bool visible);
    void layout(int x, int y, int width, int height);
    void reset(const std::string& systemLabel, const std::string& jobJson);
    void loadJobJson(const std::string& jobJson);
    std::string jobJson() const;
    std::string jobName() const;
    std::optional<SolveSystemAddress> selectedSystem() const;
    std::optional<SolveNodeAddress> selectedNode() const;
    std::optional<SolveSystemAddress> firstSystem() const;
    const SolveSystemDraft* system(const SolveSystemAddress& address) const;
    const SolveBatchDraft* batch(std::size_t index) const;
    void updateSystem(const SolveSystemAddress& address, std::string label,
                      std::string jobJson, std::size_t trajectories,
                      std::size_t clusters);
    SolveSystemAddress appendBatch(std::string batchLabel,
        std::string systemLabel, std::string jobJson,
        std::size_t trajectories, std::size_t clusters);
    bool command(int id);
    HWND treeHandle() const { return tree_; }
    bool rebuilding() const { return rebuilding_; }

private:
    HWND parent_ = nullptr;
    HWND jobCaption_ = nullptr;
    HWND jobName_ = nullptr;
    HWND tree_ = nullptr;
    HWND campaign_ = nullptr;
    HWND addBatch_ = nullptr;
    HWND addSystem_ = nullptr;
    HWND addCluster_ = nullptr;
    HWND addTrajectory_ = nullptr;
    HWND remove_ = nullptr;
    SolveJobModel model_;
    std::vector<SolveNodeAddress> visibleNodes_;
    bool rebuilding_ = false;
    void rebuild(const std::optional<SolveNodeAddress>& selection = std::nullopt);
    std::size_t selectedBatch() const;
};

}  // namespace dsl::compute
