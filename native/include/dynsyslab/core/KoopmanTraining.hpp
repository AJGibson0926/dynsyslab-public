#pragma once

#include <cstddef>
#include <vector>

namespace dynsyslab::core {

struct KoopmanTrajectoryEnsemble {
    std::vector<double> time;    // one common trajectory time vector
    std::vector<double> states;  // trajectory x sample x state
    std::size_t samples{};
    std::size_t trajectories{};
    std::size_t stateDimension{};
};

}  // namespace dynsyslab::core
