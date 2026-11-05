#include "dynsyslab_compute/SolveJobPane.hpp"

#include "dynsyslab_compute/PanelCommon.hpp"
#include "dynsyslab_viewer/Json.hpp"
#include "dynsyslab_viewer/Model.hpp"

#include <commctrl.h>

#include <algorithm>
#include <stdexcept>

namespace dsl::compute {
namespace {

HTREEITEM insertNode(HWND tree, HTREEITEM parent, const std::wstring& label,
                     LPARAM data, bool bold = false) {
    TVINSERTSTRUCTW item{};
    item.hParent = parent;
    item.hInsertAfter = TVI_LAST;
    item.item.mask = TVIF_TEXT | TVIF_PARAM | TVIF_STATE;
    item.item.pszText = const_cast<wchar_t*>(label.c_str());
    item.item.lParam = data;
    item.item.stateMask = TVIS_BOLD;
    item.item.state = bold ? TVIS_BOLD : 0;
    return TreeView_InsertItem(tree, &item);
}

bool sameNode(const SolveNodeAddress& left, const SolveNodeAddress& right) {
    return left.kind == right.kind && left.batch == right.batch
        && left.system == right.system && left.cluster == right.cluster
        && left.trajectory == right.trajectory;
}

std::size_t count(const dsl::json::Value& value, const char* key) {
    return static_cast<std::size_t>(std::max(1.0, value.at(key).number(1.0)));
}

}  // namespace

bool SolveJobPane::create(HWND parent, const dsl::theme::Fonts& fonts,
                          int campaignCommandId) {
    parent_ = parent;
    jobCaption_ = makeLabel(parent, L"Job name", fonts.caption);
    jobName_ = makeEdit(parent, kSolveJobNameId, L"New Job", fonts.body);
    tree_ = CreateWindowExW(0, WC_TREEVIEWW, L"", WS_CHILD | WS_VISIBLE
        | TVS_HASBUTTONS | TVS_HASLINES | TVS_LINESATROOT | TVS_SHOWSELALWAYS
        | TVS_FULLROWSELECT, 0, 0, 0, 0, parent,
        reinterpret_cast<HMENU>(static_cast<INT_PTR>(kSolveTreeId)),
        GetModuleHandleW(nullptr), nullptr);
    campaign_ = makeButton(parent, campaignCommandId, L"Build campaign…", fonts.body);
    addBatch_ = makeButton(parent, kAddBatchId, L"+ Batch", fonts.body);
    addSystem_ = makeButton(parent, kAddSystemId, L"+ System", fonts.body);
    addCluster_ = makeButton(parent, kAddClusterId, L"+ Cluster", fonts.body);
    addTrajectory_ = makeButton(parent, kAddTrajectoryId, L"+ Trajectory", fonts.body);
    remove_ = makeButton(parent, kRemoveSolveNodeId, L"Remove selected", fonts.body);
    if (!tree_) return false;
    dsl::theme::styleControl(tree_, fonts.body);
    TreeView_SetBkColor(tree_, dsl::theme::panel);
    TreeView_SetTextColor(tree_, dsl::theme::text);
    TreeView_SetLineColor(tree_, dsl::theme::border);
    return true;
}

void SolveJobPane::show(bool visible) {
    const int command = visible ? SW_SHOW : SW_HIDE;
    for (HWND window : {jobCaption_, jobName_, tree_, campaign_, addBatch_, addSystem_,
                        addCluster_, addTrajectory_, remove_})
        ShowWindow(window, command);
}

void SolveJobPane::layout(int x, int y, int width, int height) {
    const int gap = scalePixels(parent_, 6);
    const int captionHeight = scalePixels(parent_, 15);
    const int editHeight = scalePixels(parent_, 27);
    const int buttonHeight = scalePixels(parent_, 28);
    MoveWindow(jobCaption_, x, y, width, captionHeight, TRUE);
    MoveWindow(jobName_, x, y + captionHeight, width, editHeight, TRUE);
    const int treeY = y + captionHeight + editHeight + gap;
    const int controlsHeight = 4 * buttonHeight + 3 * gap;
    const int treeHeight = std::max(scalePixels(parent_, 100),
        height - (treeY - y) - controlsHeight - gap);
    MoveWindow(tree_, x, treeY, width, treeHeight, TRUE);
    const int half = (width - gap) / 2;
    int row = treeY + treeHeight + gap;
    MoveWindow(campaign_, x, row, width, buttonHeight, TRUE);
    row += buttonHeight + gap;
    MoveWindow(addBatch_, x, row, half, buttonHeight, TRUE);
    MoveWindow(addSystem_, x + half + gap, row, width - half - gap, buttonHeight, TRUE);
    row += buttonHeight + gap;
    MoveWindow(addCluster_, x, row, half, buttonHeight, TRUE);
    MoveWindow(addTrajectory_, x + half + gap, row, width - half - gap, buttonHeight, TRUE);
    row += buttonHeight + gap;
    MoveWindow(remove_, x, row, width, buttonHeight, TRUE);
}

void SolveJobPane::reset(const std::string& systemLabel, const std::string& jobJson) {
    model_.reset(systemLabel, jobJson);
    SetWindowTextW(jobName_, L"New Job");
    rebuild(SolveNodeAddress{SolveNodeKind::system, 0, 0, 0, 0});
}

std::string SolveJobPane::jobName() const { return utf8Text(jobName_); }

const SolveSystemDraft* SolveJobPane::system(const SolveSystemAddress& address) const {
    return model_.system(address);
}

const SolveBatchDraft* SolveJobPane::batch(std::size_t index) const {
    return model_.batch(index);
}

std::optional<SolveSystemAddress> SolveJobPane::firstSystem() const {
    return model_.firstSystem();
}

void SolveJobPane::updateSystem(const SolveSystemAddress& address, std::string label,
                                std::string jobJson, std::size_t trajectories,
                                std::size_t clusters) {
    auto* draft = model_.system(address);
    if (!draft) return;
    draft->label = std::move(label);
    draft->jobJson = std::move(jobJson);
    draft->trajectories = std::max<std::size_t>(1, trajectories);
    draft->clusters = std::clamp<std::size_t>(clusters, 1, draft->trajectories);
    rebuild(selectedNode());
}

SolveSystemAddress SolveJobPane::appendBatch(std::string batchLabel,
    std::string systemLabel, std::string jobJson,
    std::size_t trajectories, std::size_t clusters) {
    const std::size_t batch = model_.addBatch();
    model_.batches()[batch].label = std::move(batchLabel);
    const auto address = model_.addSystem(batch, std::move(systemLabel));
    auto* draft = model_.system(address);
    draft->jobJson = std::move(jobJson);
    draft->trajectories = std::max<std::size_t>(1, trajectories);
    draft->clusters = std::clamp<std::size_t>(clusters, 1, draft->trajectories);
    rebuild(SolveNodeAddress{SolveNodeKind::system, address.batch,
        address.system, 0, 0});
    return address;
}

void SolveJobPane::rebuild(const std::optional<SolveNodeAddress>& selection) {
    rebuilding_ = true;
    TreeView_DeleteAllItems(tree_);
    visibleNodes_.clear();
    HTREEITEM selectedItem = nullptr;
    auto add = [&](HTREEITEM parent, const std::wstring& label,
                   const SolveNodeAddress& address, bool bold) {
        visibleNodes_.push_back(address);
        HTREEITEM item = insertNode(tree_, parent, label,
            static_cast<LPARAM>(visibleNodes_.size()), bold);
        if (selection && sameNode(*selection, address)) selectedItem = item;
        return item;
    };
    const auto& batches = model_.batches();
    for (std::size_t batchIndex = 0; batchIndex < batches.size(); ++batchIndex) {
        const auto& batch = batches[batchIndex];
        const std::wstring batchLabel = dsl::wide(batch.label) + L"  ·  "
            + std::to_wstring(batch.systems.size()) + L" systems";
        HTREEITEM batchItem = add(TVI_ROOT, batchLabel,
            {SolveNodeKind::batch, batchIndex, 0, 0, 0}, true);
        for (std::size_t systemIndex = 0; systemIndex < batch.systems.size(); ++systemIndex) {
            const auto& system = batch.systems[systemIndex];
            std::wstring systemLabel = L"[SYSTEM] " + dsl::wide(system.label) + L"  ·  "
                + std::to_wstring(system.clusters) + L" clusters  ·  "
                + std::to_wstring(system.trajectories) + L" trajectories";
            HTREEITEM systemItem = add(batchItem, systemLabel,
                {SolveNodeKind::system, batchIndex, systemIndex, 0, 0}, true);
            for (std::size_t cluster = 0; cluster < system.clusters; ++cluster) {
                std::size_t members = 0;
                for (std::size_t trajectory = 0; trajectory < system.trajectories; ++trajectory)
                    if (assignedCluster(trajectory, system.trajectories,
                                        system.clusters) == cluster) ++members;
                HTREEITEM clusterItem = add(systemItem,
                    L"Cluster " + std::to_wstring(cluster + 1) + L"  ·  "
                        + std::to_wstring(members) + L" trajectories",
                    {SolveNodeKind::cluster, batchIndex, systemIndex, cluster, 0}, false);
                for (std::size_t trajectory = 0; trajectory < system.trajectories; ++trajectory) {
                    const std::size_t owner = assignedCluster(
                        trajectory, system.trajectories, system.clusters);
                    if (owner == cluster)
                        add(clusterItem, L"Trajectory " + std::to_wstring(trajectory + 1),
                            {SolveNodeKind::trajectory, batchIndex, systemIndex,
                             cluster, trajectory}, false);
                }
                TreeView_Expand(tree_, clusterItem, TVE_EXPAND);
            }
            TreeView_Expand(tree_, systemItem, TVE_EXPAND);
        }
        TreeView_Expand(tree_, batchItem, TVE_EXPAND);
    }
    if (!selectedItem) {
        if (const auto first = model_.firstSystem()) {
            for (std::size_t index = 0; index < visibleNodes_.size(); ++index)
                if (visibleNodes_[index].kind == SolveNodeKind::system
                    && visibleNodes_[index].batch == first->batch
                    && visibleNodes_[index].system == first->system) {
                    TVITEMW sought{}; sought.mask = TVIF_PARAM;
                    sought.lParam = static_cast<LPARAM>(index + 1);
                    HTREEITEM cursor = TreeView_GetRoot(tree_);
                    while (cursor && !selectedItem) {
                        TVITEMW current{}; current.mask = TVIF_PARAM; current.hItem = cursor;
                        TreeView_GetItem(tree_, &current);
                        if (current.lParam == sought.lParam) selectedItem = cursor;
                        cursor = TreeView_GetNextVisible(tree_, cursor);
                    }
                    break;
                }
        }
        if (!selectedItem) selectedItem = TreeView_GetRoot(tree_);
    }
    if (selectedItem) TreeView_SelectItem(tree_, selectedItem);
    rebuilding_ = false;
}

std::optional<SolveNodeAddress> SolveJobPane::selectedNode() const {
    const HTREEITEM item = TreeView_GetSelection(tree_);
    if (!item) return std::nullopt;
    TVITEMW value{}; value.mask = TVIF_PARAM; value.hItem = item;
    if (!TreeView_GetItem(tree_, &value) || value.lParam <= 0
        || static_cast<std::size_t>(value.lParam) > visibleNodes_.size()) return std::nullopt;
    return visibleNodes_[static_cast<std::size_t>(value.lParam - 1)];
}

std::optional<SolveSystemAddress> SolveJobPane::selectedSystem() const {
    const auto node = selectedNode();
    return node ? model_.owningSystem(*node) : std::nullopt;
}

std::size_t SolveJobPane::selectedBatch() const {
    const auto node = selectedNode();
    return node ? node->batch : (model_.batches().empty() ? 0 : model_.batches().size() - 1);
}

bool SolveJobPane::command(int id) {
    if (id == kAddBatchId) {
        const std::size_t batch = model_.addBatch();
        rebuild(SolveNodeAddress{SolveNodeKind::batch, batch, 0, 0, 0});
        return true;
    }
    if (id == kAddSystemId) {
        const auto address = model_.addSystem(selectedBatch());
        rebuild(SolveNodeAddress{SolveNodeKind::system, address.batch, address.system, 0, 0});
        return true;
    }
    const auto systemAddress = selectedSystem();
    if (id == kAddClusterId && systemAddress && model_.addCluster(*systemAddress)) {
        const auto* draft = model_.system(*systemAddress);
        rebuild(SolveNodeAddress{SolveNodeKind::cluster, systemAddress->batch,
            systemAddress->system, draft->clusters - 1, 0});
        return true;
    }
    if (id == kAddTrajectoryId && systemAddress && model_.addTrajectory(*systemAddress)) {
        const auto* draft = model_.system(*systemAddress);
        const std::size_t trajectory = draft->trajectories - 1;
        const std::size_t cluster = assignedCluster(
            trajectory, draft->trajectories, draft->clusters);
        rebuild(SolveNodeAddress{SolveNodeKind::trajectory, systemAddress->batch,
            systemAddress->system, cluster, trajectory});
        return true;
    }
    if (id == kRemoveSolveNodeId) {
        const auto node = selectedNode();
        if (node && model_.remove(*node)) { rebuild(); return true; }
    }
    return false;
}

void SolveJobPane::loadJobJson(const std::string& jobJson) {
    const auto root = dsl::json::parse(jobJson);
    const std::string kind = root.at("kind").string();
    model_.batches().clear();
    SetWindowTextW(jobName_, dsl::wide(root.at("label").string("New Job")).c_str());
    if (kind == "solve") {
        SolveBatchDraft batch{"Batch 1", {}};
        const auto& initial = root.at("initial_conditions");
        batch.systems.push_back({root.at("ui_system_label").string(
            root.at("system").at("name").string("System")), jobJson,
            count(initial, "count"), count(initial, "clusters")});
        model_.batches().push_back(std::move(batch));
    } else if (kind == "solve_set") {
        for (const auto& batchValue : root.at("batches").array()) {
            SolveBatchDraft batch;
            batch.label = batchValue.at("label").string(
                "Batch " + std::to_string(model_.batches().size() + 1));
            for (const auto& child : batchValue.at("systems").array()) {
                const auto& initial = child.at("initial_conditions");
                batch.systems.push_back({child.at("ui_system_label").string(
                    child.at("system").at("name").string("System")),
                    dsl::json::dump(child), count(initial, "count"),
                    count(initial, "clusters")});
            }
            model_.batches().push_back(std::move(batch));
        }
    } else throw std::invalid_argument("Solver tree requires a solve or solve_set job");
    if (model_.batches().empty()) model_.batches().push_back({"Batch 1", {}});
    rebuild();
}

std::string SolveJobPane::jobJson() const {
    dsl::json::Value::Object root;
    root["schema"] = "dynsyslab.compute-job/v1";
    root["kind"] = "solve_set";
    root["label"] = jobName().empty() ? "New Job" : jobName();
    dsl::json::Value::Array batches;
    for (const auto& batch : model_.batches()) {
        dsl::json::Value::Object encoded;
        encoded["label"] = batch.label;
        dsl::json::Value::Array systems;
        for (const auto& draft : batch.systems) {
            if (draft.jobJson.empty())
                throw std::invalid_argument("Every proposed System must be configured before planning");
            auto child = dsl::json::parse(draft.jobJson);
            child.at("schema") = "dynsyslab.compute-job/v1";
            child.at("kind") = "solve";
            child.at("label") = jobName() + " / " + batch.label + " / " + draft.label;
            child.at("ui_system_label") = draft.label;
            child.at("initial_conditions").at("count") = static_cast<double>(draft.trajectories);
            child.at("initial_conditions").at("clusters") = static_cast<double>(draft.clusters);
            systems.push_back(std::move(child));
        }
        encoded["systems"] = dsl::json::Value(std::move(systems));
        batches.emplace_back(std::move(encoded));
    }
    root["batches"] = dsl::json::Value(std::move(batches));
    return dsl::json::dump(dsl::json::Value(std::move(root)));
}

}  // namespace dsl::compute
