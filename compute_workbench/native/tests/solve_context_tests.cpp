#include "dynsyslab_compute/SolveContextPanel.hpp"

#include <iostream>
#include <stdexcept>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
}

int main() {
    try {
        const std::string duffing = R"({
          "kind":"solve","system":{"name":"duffing","preset":"double_well_reference","parameters":{}},
          "time":{"coordinate_system":"automatic"},
          "initial_conditions":{"mode":"explicit","count":4,"clusters":2,"seed":0,
            "values":[[-1,0],[-0.5,0],[0.5,0],[1,0]],"ranges":[]},
          "controller":{"kind":"none"}})";
        const auto phase = dsl::compute::clusterInitialConditionPayload(duffing, 1);
        require(phase->view == "phase" && phase->figures.size() == 1,
                "2D phase payload");
        require(phase->figures[0].square && phase->figures[0].traces[0].x.size() == 2,
                "square clustered 2D IC plot");

        const std::string lorenz = R"({
          "kind":"solve","system":{"name":"lorenz","preset":null,"parameters":{}},
          "time":{"coordinate_system":"automatic"},
          "initial_conditions":{"mode":"explicit","count":2,"clusters":1,"seed":0,
            "values":[[1,2,3],[4,5,6]],"ranges":[]},
          "controller":{"kind":"none"}})";
        const auto three = dsl::compute::clusterInitialConditionPayload(lorenz, 0);
        require(three->view == "three_d" && three->figures[0].traces[0].z.size() == 2,
                "3D clustered IC payload");
        std::cout << "dynsyslab_solve_context_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_solve_context_tests: FAIL: " << error.what() << '\n';
        return 1;
    }
}
