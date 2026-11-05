#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/workflow/Plan.hpp>
#include <dynsyslab/workflow/SystemResolution.hpp>

#include <iostream>
#include <stdexcept>
#include <string>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

dynsyslab::core::json::Value solveJob() {
    return dynsyslab::core::json::parse(R"({
      "schema":"dynsyslab.compute-job/v1","kind":"solve","label":"plan-test",
      "system":{"name":"duffing","preset":"double_well_reference","parameters":{}},
      "time":{"start":0,"stop":1,"output_step":0.01,"coordinate_system":"automatic"},
      "initial_conditions":{"mode":"default","count":2,"clusters":2},
      "controller":{"kind":"lqr","settings":{},"target":{"kind":"constant","state":[0,0]}},
      "campaign":{"enabled":true,"mode":"cartesian","axes":[
        {"path":"controller.settings.R.0.0","values":[1,2,3]}]},
      "solver":{"method":"RK45"},"archive":{"enabled":true}
    })");
}
}

int main() {
    try {
        auto job = solveJob();
        const auto plan = dynsyslab::workflow::planJob(job);
        require(plan.at("schema").string() == "dynsyslab.workflow-plan/v1", "schema");
        require(plan.at("system_count").number() == 3.0, "campaign count");
        require(plan.at("trajectory_count").number() == 6.0, "trajectory count");
        require(plan.at("cluster_count").number() == 6.0, "cluster count");
        require(plan.at("output_samples_per_trajectory").number() == 101.0, "samples");
        require(plan.at("has_control_target").boolean(), "LQR target plan metadata");
        const double feedbackBytes = plan.at("estimated_peak_solution_bytes").number();
        dynsyslab::core::SystemConfiguration delay;delay.identifier="mackey_glass";
        require(dynsyslab::workflow::resolveSystemDefinition(delay).configuration.parameters.at("delay")==17.0,
            "DDE configuration was not resolved");
        dynsyslab::core::SystemConfiguration pde;pde.identifier="heat_equation_1d";
        require(dynsyslab::workflow::resolveSystemDefinition(pde).configuration.parameters.at("points")==64.0,
            "PDE configuration was not resolved");
        dynsyslab::core::SystemConfiguration lattice;lattice.identifier="ising";
        require(dynsyslab::workflow::resolveSystemDefinition(lattice).configuration.parameters.at("size")==32.0,
            "lattice configuration was not resolved");
        job.at("controller").at("kind") = "waveform";
        job.at("controller").at("target") = dynsyslab::core::json::Value{};
        const auto waveformPlan = dynsyslab::workflow::planJob(job);
        require(!waveformPlan.at("has_control_target").boolean(),
                "prescribed waveform must be planned without target storage");
        require(waveformPlan.at("estimated_peak_solution_bytes").number() * 3.0
                == feedbackBytes * 2.0,
                "waveform storage estimate must omit state-target history");
        job.at("controller").at("kind") = "unknown";
        bool rejected = false;
        try { (void)dynsyslab::workflow::planJob(job); }
        catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "unknown controller accepted");
        std::cout << "dynsyslab_workflow_plan_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_workflow_plan_tests: FAIL: " << error.what() << '\n';
        return 1;
    }
}
