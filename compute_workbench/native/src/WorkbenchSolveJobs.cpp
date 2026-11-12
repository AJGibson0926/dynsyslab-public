#include "dynsyslab_compute/Workbench.hpp"

#include "dynsyslab_compute/PanelCommon.hpp"

#include <limits>

namespace dsl::compute {

void Workbench::storeActiveSolve() {
    if (!activeSolveSystem_) return;
    const auto* draft = solveJobs_.system(*activeSolveSystem_);
    if (!draft) { activeSolveSystem_.reset(); return; }
    const std::string label = solveJobs_.jobName().empty()
        ? "New Job" : solveJobs_.jobName();
    solveJobs_.updateSystem(*activeSolveSystem_, solver_.systemLabel(),
        solver_.jobJson(label), solver_.trajectoryCount(), solver_.clusterCount());
}

void Workbench::selectSolveNode() {
    if (!solverActive_) return;
    storeActiveSolve();
    const auto node = solveJobs_.selectedNode();
    if (!node) return;
    const auto selected = solveJobs_.selectedSystem();
    if (selected && (!activeSolveSystem_ || *activeSolveSystem_ != *selected)) {
        const auto* draft = solveJobs_.system(*selected);
        if (!draft) return;
        if (draft->jobJson.empty()) solver_.resetToSystem(
            std::numeric_limits<std::size_t>::max(), false);
        else solver_.loadJobJson(draft->jobJson, false);
        solver_.setHierarchyCounts(draft->trajectories, draft->clusters);
        activeSolveSystem_ = selected;
        // Persist hierarchy counts and normalize any per-cluster IC specs before
        // a newly added Cluster is previewed.
        storeActiveSolve();
    }
    if (node->kind == SolveNodeKind::system) {
        solveContext_.show(false); solver_.showEditor(true);
        SetWindowTextW(status_, (L"Editing proposed System: "
            + dsl::wide(solver_.systemLabel())).c_str());
    } else {
        solver_.showEditor(false);
        const auto* batch = solveJobs_.batch(node->batch);
        const auto* draft = selected ? solveJobs_.system(*selected) : nullptr;
        solveContext_.setSelection(*node,batch,draft,solveJobs_.jobName());
        solveContext_.show(true);
        const wchar_t* level=node->kind==SolveNodeKind::batch?L"Batch"
            :node->kind==SolveNodeKind::cluster?L"Cluster":L"Trajectory";
        SetWindowTextW(status_,(std::wstring(L"Inspecting proposed ")+level).c_str());
    }
}

void Workbench::editSolveTree(int id) {
    if (!solverActive_) return;
    storeActiveSolve();
    if (!solveJobs_.command(id)) {
        SetWindowTextW(status_, L"Select the appropriate parent node for that operation");
        return;
    }
    activeSolveSystem_.reset();
    selectSolveNode();
}

void Workbench::applyClusterInitialConditions() {
    if (!solverActive_) return;
    const auto address = solveJobs_.selectedSystem();
    const auto node = solveJobs_.selectedNode();
    if (!address || !node || node->kind != SolveNodeKind::cluster) return;
    const auto* current = solveJobs_.system(*address);
    if (!current) return;
    const std::string label = current->label;
    const std::size_t trajectories = current->trajectories;
    const std::size_t clusters = current->clusters;
    try {
        const auto updated = solveContext_.applyClusterConfiguration();
        if (!updated) return;
        solver_.loadJobJson(*updated, false);
        solver_.setHierarchyCounts(trajectories, clusters);
        solveJobs_.updateSystem(*address, label, *updated, trajectories, clusters);
        activeSolveSystem_ = address;
        selectSolveNode();
        SetWindowTextW(status_, L"Cluster IC generator updated in the proposed job");
    } catch (const std::exception& error) {
        SetWindowTextW(output_, dsl::wide(error.what()).c_str());
        SetWindowTextW(status_, L"Cluster IC configuration is invalid");
        MessageBoxW(window_, dsl::wide(error.what()).c_str(),
            L"Initial-condition generator", MB_OK | MB_ICONERROR);
    }
}

std::string Workbench::activeJobJson() {
    if (!solverActive_) return post_.jobJson();
    storeActiveSolve();
    return solver_.applyExecutionContract(solveJobs_.jobJson());
}

void Workbench::openKromTrainingDomainPlanner() {
    if (taskActive_) return;
    if (!solver_.prepareKromTrainingDomainExpansion(post_.source())) return;
    const std::string label = solveJobs_.jobName().empty()
        ? "New Job" : solveJobs_.jobName();
    activeSolveSystem_ = solveJobs_.appendBatch(
        "KROM training-domain expansion", solver_.systemLabel(),
        solver_.jobJson(label), solver_.trajectoryCount(),
        solver_.clusterCount());
    setTab(true);
    SetWindowTextW(status_,
        L"Prepared a new free-training solve for KROM domain expansion");
}

}  // namespace dsl::compute
