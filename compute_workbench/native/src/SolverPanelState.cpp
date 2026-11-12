#include "dynsyslab_compute/SolverPanel.hpp"

#include "dynsyslab_compute/PanelCommon.hpp"
#include "dynsyslab_viewer/Json.hpp"

#include <algorithm>
#include <stdexcept>

namespace dsl::compute {
namespace {

int comboChoice(HWND combo) {
    return static_cast<int>(SendMessageW(combo, CB_GETCURSEL, 0, 0));
}

void selectComboText(HWND combo, const std::string& text) {
    const auto value = dsl::wide(text);
    const LRESULT index = SendMessageW(combo, CB_FINDSTRINGEXACT,
        static_cast<WPARAM>(-1), reinterpret_cast<LPARAM>(value.c_str()));
    if (index != CB_ERR) SendMessageW(combo, CB_SETCURSEL, index, 0);
}

void setJson(HWND window, const dsl::json::Value& value) {
    SetWindowTextW(window, dsl::wide(dsl::json::dump(value)).c_str());
}

std::size_t positiveCount(HWND window) {
    try { return std::max<std::size_t>(1, std::stoull(utf8Text(window))); }
    catch (...) { return 1; }
}

}  // namespace

std::size_t SolverPanel::systemIndex() const {
    const int index = comboChoice(system_);
    return index < 0 ? 0 : static_cast<std::size_t>(index);
}

std::string SolverPanel::systemLabel() const {
    const std::size_t index = systemIndex();
    return index < registry_->systems.size() ? registry_->systems[index].label : "System";
}

std::size_t SolverPanel::trajectoryCount() const { return positiveCount(count_); }
std::size_t SolverPanel::clusterCount() const {
    return std::min(trajectoryCount(), positiveCount(clusters_));
}

void SolverPanel::setHierarchyCounts(std::size_t trajectories, std::size_t clusters) {
    trajectories = std::max<std::size_t>(1, trajectories);
    clusters = std::clamp<std::size_t>(clusters, 1, trajectories);
    SetWindowTextW(count_, std::to_wstring(trajectories).c_str());
    SetWindowTextW(clusters_, std::to_wstring(clusters).c_str());
}

void SolverPanel::resetToSystem(std::size_t index, bool includeExecution) {
    if (registry_->systems.empty()) throw std::runtime_error("No dynamical systems are registered");
    if (index == std::numeric_limits<std::size_t>::max()) {
        index = 0;
        for (std::size_t candidate = 0; candidate < registry_->systems.size(); ++candidate)
            if (registry_->systems[candidate].name == "rpe") { index = candidate; break; }
    }
    index = std::min(index, registry_->systems.size() - 1);
    SendMessageW(system_, CB_SETCURSEL, static_cast<WPARAM>(index), 0);
    SetWindowTextW(parameters_, L"{}");
    updatePresets(includeExecution);
    SendMessageW(preset_, CB_SETCURSEL, 0, 0);
    SetWindowTextW(start_, L"0"); SetWindowTextW(stop_, L"10");
    SetWindowTextW(step_, L"0.01");
    selectComboText(icMode_, "default"); setHierarchyCounts(1, 1);
    SetWindowTextW(seed_, L"0"); SetWindowTextW(icData_, L"[]");
    clusterSpecs_ = "[]";
    selectComboText(controller_, "none"); SetWindowTextW(controllerSettings_, L"{}");
    controlPageActive_ = false;
    updateSystemChoices(); updateControllerControls();
    SendMessageW(targetKind_, CB_SETCURSEL, 0, 0);
    SendMessageW(target_, CB_SETCURSEL, 0, 0);
    SendMessageW(targetSource_, CB_SETCURSEL, 0, 0);
    std::wstring target;
    for (int component = 0; component < std::max(1, registry_->systems[index].stateDimension);
         ++component) target += component ? L", 0" : L"0";
    SetWindowTextW(targetState_, target.c_str());
    SetWindowTextW(q_, L"1.0"); SetWindowTextW(r_, L"1.0");
    SetWindowTextW(lqrA_, L""); SetWindowTextW(lqrB_, L"");
    SetWindowTextW(targetFrequencies_, L"1");
    SetWindowTextW(targetAmplitudes_, L"0.1");
    SetWindowTextW(targetPhases_, L"0");
    if (includeExecution)
        campaign_ = "{\"enabled\":false,\"mode\":\"cartesian\",\"axes\":[]}";
    selectComboText(representation_, "exact");
    selectComboText(controlEvaluation_, "continuous RHS");
    SetWindowTextW(controlSampleTime_, L"");
    SetWindowTextW(horizon_, L"20");
    SetWindowTextW(controlLimit_, L""); clearTrainingPreflight();
    updateParameterSummaries(); updatePanelPageVisibility();
}

void SolverPanel::loadJobJson(const std::string& jobJson, bool includeExecution) {
    const auto root = dsl::json::parse(jobJson);
    if (root.at("kind").string() != "solve")
        throw std::invalid_argument("Solver form requires a solve job");
    const std::string systemName = root.at("system").at("name").string();
    std::size_t index = 0;
    for (; index < registry_->systems.size(); ++index)
        if (registry_->systems[index].name == systemName) break;
    if (index == registry_->systems.size())
        throw std::invalid_argument("Saved solve references an unavailable dynamical system");
    resetToSystem(index, includeExecution);
    selectComboText(preset_, root.at("system").at("preset").string());
    setJson(parameters_, root.at("system").at("parameters"));
    const auto& time = root.at("time");
    setJson(start_, time.at("start")); setJson(stop_, time.at("stop"));
    setJson(step_, time.at("output_step"));
    std::string coordinates = time.at("coordinate_system").string();
    if (coordinates.empty() || coordinates == "automatic")
        coordinates = registry_->systems[index].preferredCoordinates.empty()
            ? "dimensional" : registry_->systems[index].preferredCoordinates;
    selectComboText(coordinates_, coordinates);
    updateTimeLabels();
    const auto& solver = root.at("solver");
    const std::string evaluation = solver.at("control_evaluation").string();
    selectComboText(controlEvaluation_, evaluation == "sampled"
        ? "sampled" : "continuous RHS");
    if (!solver.at("control_sample_time").isNull())
        setJson(controlSampleTime_, solver.at("control_sample_time"));
    else SetWindowTextW(controlSampleTime_, L"");
    const auto& canonicalPartition = root.at("system").at("pde_partitioning");
    const auto& legacyPartition = solver.at("parallel");
    const auto& partition = canonicalPartition.isNull()
        ? legacyPartition : canonicalPartition;
    SendMessageW(pdePartition_, BM_SETCHECK,
        partition.at("enabled").boolean(false) ? BST_CHECKED : BST_UNCHECKED, 0);
    if (!partition.at("partition_shape").isNull())
        setJson(pdePartitionShape_, partition.at("partition_shape"));
    if (!partition.at("ghost_width").isNull())
        setJson(pdeGhostWidth_, partition.at("ghost_width"));
    if (!partition.at("local_workers").isNull())
        setJson(pdeLocalWorkers_, partition.at("local_workers"));
    const auto& initial = root.at("initial_conditions");
    clusterSpecs_ = dsl::json::dump(initial.at("cluster_specs"));
    if (initial.at("cluster_specs").isNull()) clusterSpecs_ = "[]";
    const std::string mode = initial.at("mode").string("default");
    selectComboText(icMode_, mode);
    setHierarchyCounts(static_cast<std::size_t>(std::max(1.0, initial.at("count").number(1.0))),
        static_cast<std::size_t>(std::max(1.0, initial.at("clusters").number(1.0))));
    setJson(seed_, initial.at("seed"));
    setJson(icData_, mode == "explicit" ? initial.at("values")
        : (mode == "sphere" || mode == "ellipsoid")
            ? initial.at("distribution") : initial.at("ranges"));
    const auto& controller = root.at("controller");
    selectComboText(controller_, controller.at("kind").string("none"));
    controlPageActive_ = false;
    updateControllerControls();
    const auto& settings = controller.at("settings");
    setJson(controllerSettings_, settings);
    if (!settings.at("Q").isNull()) setJson(q_, settings.at("Q"));
    if (!settings.at("R").isNull()) setJson(r_, settings.at("R"));
    selectComboText(representation_, settings.at("control_representation").string("exact"));
    if (!settings.at("horizon").isNull()) setJson(horizon_, settings.at("horizon"));
    if (!settings.at("control_limit").isNull()) setJson(controlLimit_, settings.at("control_limit"));
    if (!settings.at("A").isNull()) setJson(lqrA_, settings.at("A"));
    if (!settings.at("B").isNull()) setJson(lqrB_, settings.at("B"));
    const std::string krom = settings.at("krom").string();
    for (std::size_t visible = 0; visible < visibleKroms_.size(); ++visible)
        if (registry_->kroms[visibleKroms_[visible]].source == krom) {
            SendMessageW(krom_, CB_SETCURSEL, static_cast<WPARAM>(visible + 1), 0); break;
        }
    const auto& targetSpec=controller.at("target");
    const auto& target=targetSpec.at("kind").string()=="from_state"
        ?targetSpec.at("state_target"):targetSpec;
    const std::string targetKind = target.at("kind").string("constant");
    if (targetKind == "rpe_canonical") {
        SendMessageW(targetKind_, CB_SETCURSEL, 3, 0);
        selectComboText(target_, target.at("name").string());
    } else if (targetKind == "trajectory") {
        SendMessageW(targetKind_, CB_SETCURSEL, 1, 0);
        const std::string source = target.at("source").string();
        for (std::size_t visible = 0; visible < visibleTargetTrajectories_.size(); ++visible)
            if (registry_->trajectories[visibleTargetTrajectories_[visible]].source == source) {
                SendMessageW(targetSource_, CB_SETCURSEL,
                    static_cast<WPARAM>(visible + 1), 0); break;
            }
    } else if (targetKind == "waveform") {
        SendMessageW(targetKind_, CB_SETCURSEL, 2, 0);
        setJson(targetState_, target.at("offset"));
        setJson(targetFrequencies_, target.at("frequencies"));
        setJson(targetAmplitudes_, target.at("amplitudes"));
        setJson(targetPhases_, target.at("phases_rad"));
    } else {
        SendMessageW(targetKind_, CB_SETCURSEL, 0, 0);
        if (!target.isNull()) setJson(targetState_, target.at("state"));
    }
    if(targetSpec.at("space").string()=="observable")
        SendMessageW(targetKind_,CB_SETCURSEL,SendMessageW(targetKind_,CB_GETCURSEL,0,0)
            +(targetSpec.at("kind").string()=="from_state"?8:4),0);
    if (selectedText(controller_) == "waveform") updateWaveformControls();
    updatePanelPageVisibility();
    if (includeExecution) {
        campaign_ = dsl::json::dump(root.at("campaign"));
        execution_.load(dsl::json::dump(solver));
    }
    preflight_ = dsl::json::dump(root.at("preflight"));
    archiveDataRole_ = root.at("archive").at("data_role").string("free");
    updateParameterSummaries();
}

}  // namespace dsl::compute
