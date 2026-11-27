#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/workflow/Solve.hpp>

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
}

int main() {
    try {
        const auto job = dynsyslab::core::json::parse(R"({
          "kind":"solve","system":{"name":"duffing","preset":"double_well_reference","parameters":{}},
          "initial_conditions":{"mode":"default","count":4,"clusters":2,"seed":0,
            "values":[],"ranges":[],"cluster_specs":[
              {"mode":"sphere","seed":17,"values":[],"ranges":[],
               "distribution":{"center":[-1,0],"radius":0.2,"inner_radius":1,"outer_radius":1}},
              {"mode":"uniform","seed":23,"values":[],"ranges":[[0.5,1],[2,3]],
               "distribution":null}]},
          "controller":{"kind":"none"}})");
        const auto first = dynsyslab::workflow::resolveOdeInitialConditions(job);
        const auto second = dynsyslab::workflow::resolveOdeInitialConditions(job);
        require(first.trajectoryCount == 4 && first.stateDimension == 2,
                "initial-condition shape");
        require(first.values == second.values, "seeded generation is not deterministic");
        for (std::size_t row = 0; row < 2; ++row) {
            const double x = first.values[row * 2] + 1.0;
            const double y = first.values[row * 2 + 1];
            require(std::abs(std::sqrt(x * x + y * y) - 0.2) < 1.0e-12,
                    "spherical cluster radius");
        }
        for (std::size_t row = 2; row < 4; ++row)
            require(first.values[row * 2] >= 0.5 && first.values[row * 2] <= 1.0
                && first.values[row * 2 + 1] >= 2.0 && first.values[row * 2 + 1] <= 3.0,
                "uniform cluster bounds");
        const auto encoded = dynsyslab::workflow::odeInitialConditionSetJson(first);
        require(encoded.at("schema").string() == "dynsyslab.initial-conditions/v1"
            && encoded.at("values").array().size() == 4,
            "machine-readable IC response");
        require(dynsyslab::workflow::supportedOdeInitialConditionModes().size() == 12,
                "IC capability list");
        auto partial = job;
        partial.at("initial_conditions").at("cluster_specs") =
            dynsyslab::core::json::parse(R"([null,{"mode":"uniform","seed":5,
              "values":[],"ranges":[[0.5,1],[2,3]],"distribution":null}])");
        const auto inherited = dynsyslab::workflow::resolveOdeInitialConditions(partial);
        auto baselineJob = partial;
        baselineJob.at("initial_conditions").at("cluster_specs") =
            dynsyslab::core::json::parse("[]");
        const auto baseline = dynsyslab::workflow::resolveOdeInitialConditions(baselineJob);
        require(inherited.values[0] == baseline.values[0]
            && inherited.values[1] == baseline.values[1]
            && inherited.values[2] == baseline.values[2]
            && inherited.values[3] == baseline.values[3],
            "null cluster spec did not inherit global defaults");
        std::cout << "dynsyslab_initial_condition_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_initial_condition_tests: FAIL: " << error.what() << '\n';
        return 1;
    }
}
