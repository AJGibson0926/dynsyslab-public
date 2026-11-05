#include "dynsyslab_compute/SolveJobModel.hpp"

#include <iostream>
#include <stdexcept>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
}

int main() {
    try {
        dsl::compute::SolveJobModel model;
        model.reset("Duffing oscillator", "{\"kind\":\"solve\"}");
        auto first = model.firstSystem();
        require(first.has_value(), "default system missing");
        require(model.addCluster(*first), "cluster addition failed");
        require(model.system(*first)->clusters == 2
            && model.system(*first)->trajectories == 2,
            "cluster addition did not preserve a valid hierarchy");
        require(model.addTrajectory(*first)
            && model.system(*first)->trajectories == 3,
            "trajectory addition failed");
        const auto batch = model.addBatch();
        const auto second = model.addSystem(batch, "RPE");
        require(model.batches().size() == 2 && model.system(second),
                "batch/system addition failed");
        require(model.remove({dsl::compute::SolveNodeKind::system,
            second.batch, second.system, 0, 0}), "system removal failed");
        require(model.remove({dsl::compute::SolveNodeKind::batch,
            0, 0, 0, 0}) && model.batches().size() == 1,
            "batch removal produced an invalid model");
        std::cout << "dynsyslab_solve_job_model_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_solve_job_model_tests: FAIL: " << error.what() << '\n';
        return 1;
    }
}
