#include "dynsyslab_compute/SolverExecutionContract.hpp"

#include <dynsyslab/core/Json.hpp>

#include <iostream>
#include <stdexcept>
#include <string>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
}

int main() {
    try {
        const std::string job = R"({
            "kind":"solve_set","batches":[{"systems":[
                {"time":{"start":0,"stop":10,"output_step":0.01,
                    "coordinate_system":"dimensionless"},
                    "system":{"family":"ode"},
                    "solver":{"control_evaluation":"continuous_rhs"},"campaign":{}},
                {"time":{"start":2,"stop":40,"output_step":0.2,
                    "coordinate_system":"dimensional"},
                    "system":{"family":"pde"},
                    "solver":{"control_evaluation":"sampled",
                        "control_sample_time":0.02,
                        "parallel":{"enabled":true,"partition_shape":[4],
                            "ghost_width":1,"local_workers":2}},"campaign":{}}
            ]}]})";
        const auto result = dynsyslab::core::json::parse(
            dsl::compute::applySolverExecutionContract(job,
                R"({"method":"BDF","rtol":1e-9})",
                R"({"enabled":true,"mode":"cartesian","axes":[]})"));
        const auto& systems = result.at("batches").array().front()
            .at("systems").array();
        require(systems[0].at("time").at("stop").number() == 10.0,
            "first System time contract was overwritten");
        require(systems[1].at("time").at("start").number() == 2.0,
            "second System time contract was overwritten");
        require(systems[0].at("time").at("coordinate_system").string()
                == "dimensionless"
            && systems[1].at("time").at("coordinate_system").string()
                == "dimensional",
            "System coordinate contracts were overwritten");
        require(systems[0].at("solver").at("control_evaluation").string()
                == "continuous_rhs"
            && systems[1].at("solver").at("control_evaluation").string()
                == "sampled"
            && systems[1].at("solver").at("control_sample_time").number()
                == 0.02,
            "per-System Controller evaluation contract was overwritten");
        require(systems[1].at("system").at("pde_partitioning")
                .at("partition_shape").array().front().number() == 4.0
            && systems[1].at("solver").at("parallel").isNull(),
            "legacy PDE partitioning was not migrated to its System");
        for (const auto& system : systems) {
            require(system.at("solver").at("method").string() == "BDF",
                "job-wide solver settings were not applied");
            require(system.at("campaign").at("enabled").boolean(false),
                "job-wide campaign settings were not applied");
        }
        std::cout << "solver execution contract tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n'; return 1;
    }
}
