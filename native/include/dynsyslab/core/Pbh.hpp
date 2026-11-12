#pragma once
#include <complex>
#include <cstddef>
#include <vector>
namespace dynsyslab::core {
struct PbhResult {
    std::vector<std::complex<double>> eigenvalues;
    std::vector<double> relativeMargins; // sigma_min / sigma_max of each PBH matrix
    bool fullRank{true}; // observable / controllable
    bool stableHiddenModes{true}; // detectable / stabilizable
    bool discreteTime{};
};
[[nodiscard]] PbhResult pbhObservability(const std::vector<double>& A,
    const std::vector<double>& C,std::size_t states,std::size_t outputs,
    bool discreteTime=false,double relativeTolerance=1e-10);
[[nodiscard]] PbhResult pbhControllability(const std::vector<double>& A,
    const std::vector<double>& B,std::size_t states,std::size_t inputs,
    bool discreteTime=false,double relativeTolerance=1e-10);
} // namespace dynsyslab::core
