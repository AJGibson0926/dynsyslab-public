#include "dynsyslab_compute/SolveContextPanel.hpp"

#include "dynsyslab_compute/PanelCommon.hpp"
#include "dynsyslab_viewer/Json.hpp"

#include <dynsyslab/workflow/Solve.hpp>

#include <algorithm>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace dsl::compute {
namespace {

std::wstring stateText(const dynsyslab::workflow::OdeInitialConditionSet& set,
                       std::size_t trajectory) {
    std::wostringstream text;
    text << L"(" << std::setprecision(7);
    for (std::size_t axis = 0; axis < set.stateDimension; ++axis) {
        if (axis) text << L", ";
        const std::string label = axis < set.stateLabels.size()
            ? set.stateLabels[axis] : "x" + std::to_string(axis + 1);
        text << dsl::wide(label) << L" = "
             << set.values[trajectory * set.stateDimension + axis];
    }
    return text.str() + L")";
}

const dsl::json::Value& clusterSpec(const dsl::json::Value& initial,
                                    std::size_t cluster) {
    const auto& specs = initial.at("cluster_specs").array();
    return cluster < specs.size() && !specs[cluster].isNull()
        ? specs[cluster] : initial;
}

std::wstring clusterDetail(const dsl::json::Value& job,
                           const dynsyslab::workflow::OdeInitialConditionSet& set,
                           const SolveSystemDraft& draft, std::size_t cluster) {
    std::size_t members = 0;
    for (std::size_t trajectory = 0; trajectory < set.trajectoryCount; ++trajectory)
        if (assignedCluster(trajectory, set.trajectoryCount, draft.clusters) == cluster)
            ++members;
    const auto& initial = job.at("initial_conditions");
    const auto& spec = clusterSpec(initial, cluster);
    std::wostringstream detail;
    detail << L"System\r\n" << dsl::wide(draft.label)
           << L"\r\n\r\nCluster membership\r\n" << members
           << L" of " << set.trajectoryCount << L" trajectories"
           << L"\r\n\r\nIC generator\r\n"
           << dsl::wide(spec.at("mode").string("default"))
           << L"\r\n\r\nSeed\r\n" << std::setprecision(12)
           << spec.at("seed").number(0.0)
           << L"\r\n\r\nCoordinates\r\n"
           << dsl::wide(job.at("time").at("coordinate_system").string("automatic"));
    return detail.str();
}

}  // namespace

std::shared_ptr<dsl::Payload> clusterInitialConditionPayload(
    const std::string& jobJson, std::size_t cluster) {
    const auto job = dsl::json::parse(jobJson);
    const auto set = dynsyslab::workflow::resolveOdeInitialConditions(job);
    const std::size_t clusters = static_cast<std::size_t>(std::max(
        1.0, job.at("initial_conditions").at("clusters").number(1.0)));
    if (cluster >= clusters) throw std::invalid_argument("selected cluster is unavailable");
    if (set.stateDimension != 2 && set.stateDimension != 3)
        throw std::invalid_argument("IC phase preview is available for 2D and 3D ODE systems");
    auto payload = std::make_shared<dsl::Payload>();
    payload->view = set.stateDimension == 2 ? "phase" : "three_d";
    payload->tileName = "Proposed initial conditions";
    dsl::Figure figure;
    figure.kind = payload->view;
    figure.title = set.stateDimension == 2
        ? "Cluster initial conditions" : "Cluster initial conditions · 3D";
    figure.xLabel = set.stateLabels.empty() ? "x1" : set.stateLabels[0];
    figure.yLabel = set.stateLabels.size() < 2 ? "x2" : set.stateLabels[1];
    figure.square = true;
    dsl::Trace points;
    points.label = "trajectory ICs";
    points.role = "initial_conditions";
    points.color = "#4FC3F7";
    points.marker = true;
    for (std::size_t trajectory = 0; trajectory < set.trajectoryCount; ++trajectory) {
        if (assignedCluster(trajectory, set.trajectoryCount, clusters) != cluster) continue;
        const std::size_t offset = trajectory * set.stateDimension;
        points.x.push_back(set.values[offset]);
        points.y.push_back(set.values[offset + 1]);
        if (set.stateDimension == 3) points.z.push_back(set.values[offset + 2]);
    }
    figure.traces.push_back(std::move(points));
    payload->figures.push_back(std::move(figure));
    return payload;
}

bool SolveContextPanel::create(HWND parent, const dsl::theme::Fonts& fonts) {
    parent_ = parent;
    kind_ = makeLabel(parent, L"Batch", fonts.caption);
    title_ = makeLabel(parent, L"Batch 1", fonts.title);
    detail_ = makeEdit(parent, 0, L"", fonts.body, true);
    hint_ = makeLabel(parent, L"", fonts.caption);
    modeLabel_ = makeLabel(parent, L"IC generation mode", fonts.caption);
    mode_ = makeCombo(parent, kClusterIcModeId, fonts.body);
    setComboItems(mode_, dynsyslab::workflow::supportedOdeInitialConditionModes());
    seedLabel_ = makeLabel(parent, L"Random seed", fonts.caption);
    seed_ = makeEdit(parent, kClusterIcSeedId, L"0", fonts.body);
    dataLabel_ = makeLabel(parent, L"Distribution data", fonts.caption);
    data_ = makeEdit(parent, kClusterIcDataId, L"[]", fonts.mono, true);
    apply_ = makeButton(parent, kClusterIcApplyId, L"Apply and regenerate preview", fonts.body);
    SendMessageW(detail_, EM_SETREADONLY, TRUE, 0);
    return kind_ && title_ && detail_ && hint_ && modeLabel_ && mode_ && seedLabel_
        && seed_ && dataLabel_ && data_ && apply_ && plot_.create(parent, fonts);
}

void SolveContextPanel::show(bool visible) {
    visible_ = visible;
    for (HWND control : {kind_, title_, detail_, hint_})
        ShowWindow(control, visible ? SW_SHOW : SW_HIDE);
    const bool clusterControls = visible && selectionKind_ == SolveNodeKind::cluster;
    for (HWND control : {modeLabel_, mode_, seedLabel_, seed_, dataLabel_, data_, apply_})
        ShowWindow(control, clusterControls ? SW_SHOW : SW_HIDE);
    ShowWindow(plot_.handle(), visible && plotVisible_ ? SW_SHOW : SW_HIDE);
}

void SolveContextPanel::layout(int x, int y, int width, int height) {
    layoutX_ = x; layoutY_ = y; layoutWidth_ = width; layoutHeight_ = height;
    const int gap = scalePixels(parent_, 12);
    const int kindHeight = scalePixels(parent_, 18);
    const int titleHeight = scalePixels(parent_, 38);
    const int top = y + kindHeight + titleHeight + gap;
    const int available = std::max(scalePixels(parent_, 120), height - (top - y));
    int textWidth = width;
    if (selectionKind_ == SolveNodeKind::cluster && plotVisible_)
        textWidth = std::max(scalePixels(parent_, 210), width * 34 / 100);
    MoveWindow(kind_, x, y, textWidth, kindHeight, TRUE);
    MoveWindow(title_, x, y + kindHeight, textWidth, titleHeight, TRUE);
    if (selectionKind_ == SolveNodeKind::cluster) {
        const int labelH = scalePixels(parent_, 15), fieldH = scalePixels(parent_, 28);
        int row = top;
        MoveWindow(modeLabel_, x, row, textWidth, labelH, TRUE);
        moveFieldControl(mode_, x, row + labelH, textWidth, fieldH);
        row += labelH + fieldH + gap;
        MoveWindow(seedLabel_, x, row, textWidth, labelH, TRUE);
        MoveWindow(seed_, x, row + labelH, textWidth, fieldH, TRUE);
        row += labelH + fieldH + gap;
        MoveWindow(dataLabel_, x, row, textWidth, labelH, TRUE);
        const int dataH = scalePixels(parent_, 105);
        MoveWindow(data_, x, row + labelH, textWidth, dataH, TRUE);
        row += labelH + dataH + gap;
        MoveWindow(apply_, x, row, textWidth, scalePixels(parent_, 32), TRUE);
        row += scalePixels(parent_, 44);
        MoveWindow(detail_, x, row, textWidth,
            std::max(scalePixels(parent_, 70), y + height - row - scalePixels(parent_, 30)), TRUE);
    } else MoveWindow(detail_, x, top, textWidth,
        available - scalePixels(parent_, 30), TRUE);
    MoveWindow(hint_, x, top + available - scalePixels(parent_, 24), textWidth,
        scalePixels(parent_, 22), TRUE);
    if (plotVisible_) MoveWindow(plot_.handle(), x + textWidth + gap, y,
        std::max(1, width - textWidth - gap), height, TRUE);
}

void SolveContextPanel::setSelection(const SolveNodeAddress& node,
                                     const SolveBatchDraft* batch,
                                     const SolveSystemDraft* system,
                                     const std::string& jobName) {
    selectionKind_ = node.kind;
    plotVisible_ = false;
    jobJson_.clear(); cluster_ = node.cluster;
    std::wstring kind, title, detail, hint;
    if (node.kind == SolveNodeKind::batch) {
        kind = L"Batch";
        title = batch ? dsl::wide(batch->label) : L"Batch";
        std::size_t clusters = 0, trajectories = 0;
        if (batch) for (const auto& child : batch->systems) {
            clusters += child.clusters; trajectories += child.trajectories;
        }
        std::wostringstream summary;
        summary << L"Job\r\n" << dsl::wide(jobName)
                << L"\r\n\r\nSystems\r\n" << (batch ? batch->systems.size() : 0)
                << L"\r\n\r\nClusters\r\n" << clusters
                << L"\r\n\r\nTrajectories\r\n" << trajectories;
        detail = summary.str();
        hint = L"Select a System to configure its dynamics.";
    } else if (!system || system->jobJson.empty()) {
        kind = L"System"; title = L"New System";
        detail = L"This proposed System has not been configured.";
        hint = L"Select the System node to configure it.";
    } else try {
        jobJson_ = system->jobJson;
        const auto job = dsl::json::parse(jobJson_);
        const auto set = dynsyslab::workflow::resolveOdeInitialConditions(job);
        if (node.kind == SolveNodeKind::cluster) {
            kind = L"Initial-condition cluster";
            title = L"Cluster " + std::to_wstring(node.cluster + 1);
            detail = clusterDetail(job, set, *system, node.cluster);
            hint = L"Axes are fitted automatically to the proposed ICs.";
            plotVisible_ = true;
            const auto& initial = job.at("initial_conditions");
            const auto& spec = clusterSpec(initial, node.cluster);
            const std::string mode = spec.at("mode").string("default");
            const auto modeText = dsl::wide(mode);
            const LRESULT modeIndex = SendMessageW(mode_, CB_FINDSTRINGEXACT,
                static_cast<WPARAM>(-1), reinterpret_cast<LPARAM>(modeText.c_str()));
            if (modeIndex != CB_ERR) SendMessageW(mode_, CB_SETCURSEL, modeIndex, 0);
            SetWindowTextW(seed_, dsl::wide(dsl::json::dump(spec.at("seed"))).c_str());
            const auto& data = mode == "explicit" ? spec.at("values")
                : (mode == "sphere" || mode == "ellipsoid")
                    ? spec.at("distribution") : spec.at("ranges");
            SetWindowTextW(data_, dsl::wide(dsl::json::dump(data)).c_str());
            updateModeHelp();
            try { plot_.setPayload(clusterInitialConditionPayload(system->jobJson, node.cluster)); }
            catch (const std::exception& error) {
                plot_.setMessage(L"IC plot unavailable", dsl::wide(error.what()));
            }
        } else {
            kind = L"Trajectory";
            title = L"Trajectory " + std::to_wstring(node.trajectory + 1);
            detail = L"System\r\n" + dsl::wide(system->label)
                + L"\r\n\r\nCluster\r\nCluster " + std::to_wstring(node.cluster + 1)
                + L"\r\n\r\nInitial condition\r\n"
                + stateText(set, node.trajectory);
            hint = L"Select the owning Cluster to inspect all of its ICs.";
        }
    } catch (const std::exception& error) {
        kind = node.kind == SolveNodeKind::cluster ? L"Initial-condition cluster" : L"Trajectory";
        title = node.kind == SolveNodeKind::cluster
            ? L"Cluster " + std::to_wstring(node.cluster + 1)
            : L"Trajectory " + std::to_wstring(node.trajectory + 1);
        detail = L"The proposed configuration could not be previewed.\r\n\r\n"
            + dsl::wide(error.what());
        hint = L"Correct the owning System or IC configuration before execution.";
        if (node.kind == SolveNodeKind::cluster) {
            plotVisible_ = true;
            plot_.setMessage(L"IC plot unavailable", dsl::wide(error.what()));
        }
    }
    SetWindowTextW(kind_, kind.c_str()); SetWindowTextW(title_, title.c_str());
    SetWindowTextW(detail_, detail.c_str()); SetWindowTextW(hint_, hint.c_str());
    if (layoutWidth_ > 0) layout(layoutX_, layoutY_, layoutWidth_, layoutHeight_);
    if (visible_) show(true);
}

std::optional<std::string> SolveContextPanel::applyClusterConfiguration() {
    if (selectionKind_ != SolveNodeKind::cluster || jobJson_.empty()) return std::nullopt;
    auto job = dsl::json::parse(jobJson_);
    auto& initial = job.at("initial_conditions");
    const std::size_t clusterCount = static_cast<std::size_t>(std::max(
        1.0, initial.at("clusters").number(1.0)));
    auto specs = initial.at("cluster_specs").array();
    while (specs.size() < clusterCount) specs.emplace_back();
    specs.resize(clusterCount);
    if (cluster_ >= specs.size()) throw std::invalid_argument("selected cluster is unavailable");
    auto& spec = specs[cluster_];
    spec = dsl::json::Value(dsl::json::Value::Object{});
    const std::string mode = selectedText(mode_);
    spec.at("mode") = mode;
    spec.at("seed") = dsl::json::parse(utf8Text(seed_));
    spec.at("values") = dsl::json::Value(dsl::json::Value::Array{});
    spec.at("ranges") = dsl::json::Value(dsl::json::Value::Array{});
    spec.at("distribution") = dsl::json::Value();
    const auto data = dsl::json::parse(utf8Text(data_).empty() ? "[]" : utf8Text(data_));
    if (mode == "explicit") spec.at("values") = data;
    else if (mode == "uniform" || mode == "normal") spec.at("ranges") = data;
    else if (mode == "sphere" || mode == "ellipsoid") spec.at("distribution") = data;
    else if (mode != "default") throw std::invalid_argument("unsupported IC generation mode");
    initial.at("cluster_specs") = dsl::json::Value(std::move(specs));
    (void)dynsyslab::workflow::resolveOdeInitialConditions(job);
    jobJson_ = dsl::json::dump(job);
    return jobJson_;
}

void SolveContextPanel::command(int id, int notification) {
    if (id == kClusterIcModeId && notification == CBN_SELCHANGE) updateModeHelp();
}

void SolveContextPanel::updateModeHelp() {
    const std::string mode = selectedText(mode_);
    const wchar_t* label = L"No parameters (system default)";
    if (mode == "explicit") label = L"Values  [[x1, x2, ...], ...]";
    else if (mode == "uniform" || mode == "normal")
        label = L"Per-state bounds  [[min, max], ...]";
    else if (mode == "sphere")
        label = L"Distribution  {center, radius, inner_radius, outer_radius}";
    else if (mode == "ellipsoid")
        label = L"Distribution  {center, axes, inner_radius, outer_radius}";
    SetWindowTextW(dataLabel_, label);
    EnableWindow(data_, mode != "default");
}

}  // namespace dsl::compute
