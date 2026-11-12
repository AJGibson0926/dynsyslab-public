#include "dynsyslab_compute/SolveJobPane.hpp"

#include "dynsyslab_viewer/Json.hpp"

#include <commctrl.h>

#include <iostream>
#include <stdexcept>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
}

int main() {
    HWND parent = nullptr;
    try {
        INITCOMMONCONTROLSEX controls{sizeof(controls), ICC_TREEVIEW_CLASSES};
        require(InitCommonControlsEx(&controls) != FALSE, "common controls unavailable");
        parent = CreateWindowExW(0, L"STATIC", L"solve-job-pane-test",
            WS_OVERLAPPEDWINDOW, 0, 0, 600, 800, nullptr, nullptr,
            GetModuleHandleW(nullptr), nullptr);
        require(parent != nullptr, "test parent creation failed");
        dsl::theme::Fonts fonts;
        dsl::compute::SolveJobPane pane;
        require(pane.create(parent, fonts, 7777), "solve job pane creation failed");
        const std::string child = R"({
          "schema":"dynsyslab.compute-job/v1","kind":"solve","label":"test",
          "system":{"name":"duffing","preset":null,"parameters":{}},
          "time":{"start":0,"stop":1,"output_step":0.1,"coordinate_system":"automatic"},
          "initial_conditions":{"mode":"default","count":1,"clusters":1,"seed":0,"values":[],"ranges":[]},
          "controller":{"kind":"none","settings":{},"target":{"kind":"constant","state":[0,0]}},
          "campaign":{"enabled":false},"preflight":{"enabled":false},
          "solver":{"method":"RK4"},"archive":{"enabled":true}}
        )";
        pane.reset("Duffing oscillator", child);
        pane.layout(0, 0, 300, 700);
        require(TreeView_GetCount(pane.treeHandle()) == 4,
                "default hierarchy should contain batch/system/cluster/trajectory");
        require(pane.command(dsl::compute::kAddTrajectoryId),
                "add trajectory command failed");
        require(TreeView_GetCount(pane.treeHandle()) == 5,
                "trajectory was not added to the visible tree");
        require(pane.command(dsl::compute::kAddSystemId),
                "add system command failed");
        const auto selected = pane.selectedSystem();
        require(selected && pane.system(*selected)->jobJson.empty(),
                "new system was not selected as an unconfigured draft");
        pane.updateSystem(*selected, "RPE", child, 2, 1);
        const auto encoded = dsl::json::parse(pane.jobJson());
        require(encoded.at("kind").string() == "solve_set",
                "pane did not encode a solve_set job");
        require(encoded.at("batches").array()[0].at("systems").array().size() == 2,
                "pane lost a configured system");
        const auto expansion = pane.appendBatch("KROM training-domain expansion",
            "Duffing oscillator", child, 12, 1);
        require(expansion.batch == 1 && expansion.system == 0
            && pane.system(expansion)->trajectories == 12,
            "KROM training-domain handoff did not append a configured Batch");
        const auto expanded = dsl::json::parse(pane.jobJson());
        require(expanded.at("batches").array().size() == 2
            && expanded.at("batches").array()[1].at("label").string()
                == "KROM training-domain expansion",
            "appended training Batch was not retained in the proposed job");
        DestroyWindow(parent);
        std::cout << "dynsyslab_solve_job_pane_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        if (parent) DestroyWindow(parent);
        std::cerr << "dynsyslab_solve_job_pane_tests: FAIL: " << error.what() << '\n';
        return 1;
    }
}
