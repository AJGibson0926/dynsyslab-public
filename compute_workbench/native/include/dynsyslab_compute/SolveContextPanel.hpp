#pragma once

#include "dynsyslab_compute/SolveJobModel.hpp"
#include "dynsyslab_viewer/PlotCanvas.hpp"
#include "dynsyslab_viewer/Theme.hpp"

#include <windows.h>

#include <memory>
#include <optional>
#include <string>

namespace dsl::compute {

constexpr int kClusterIcModeId = 5230;
constexpr int kClusterIcSeedId = 5231;
constexpr int kClusterIcDataId = 5232;
constexpr int kClusterIcApplyId = 5233;

[[nodiscard]] std::shared_ptr<dsl::Payload> clusterInitialConditionPayload(
    const std::string& jobJson, std::size_t cluster);

class SolveContextPanel {
public:
    bool create(HWND parent, const dsl::theme::Fonts& fonts);
    void show(bool visible);
    void layout(int x, int y, int width, int height);
    void setSelection(const SolveNodeAddress& node, const SolveBatchDraft* batch,
                      const SolveSystemDraft* system, const std::string& jobName);
    void command(int id, int notification);
    [[nodiscard]] std::optional<std::string> applyClusterConfiguration();

private:
    HWND parent_ = nullptr;
    HWND kind_ = nullptr;
    HWND title_ = nullptr;
    HWND detail_ = nullptr;
    HWND hint_ = nullptr;
    HWND modeLabel_ = nullptr;
    HWND mode_ = nullptr;
    HWND seedLabel_ = nullptr;
    HWND seed_ = nullptr;
    HWND dataLabel_ = nullptr;
    HWND data_ = nullptr;
    HWND apply_ = nullptr;
    dsl::PlotCanvas plot_;
    SolveNodeKind selectionKind_ = SolveNodeKind::batch;
    bool visible_ = false;
    bool plotVisible_ = false;
    int layoutX_ = 0;
    int layoutY_ = 0;
    int layoutWidth_ = 0;
    int layoutHeight_ = 0;
    std::string jobJson_;
    std::size_t cluster_ = 0;
    void updateModeHelp();
};

}  // namespace dsl::compute
