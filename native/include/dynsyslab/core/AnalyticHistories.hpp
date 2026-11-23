#pragma once
#include <dynsyslab/core/AnalyticSolution.hpp>
#include <dynsyslab/core/Observable.hpp>
namespace dynsyslab::core {
struct AnalyticHistories {
    std::size_t points{},times{},channels{};
    std::vector<std::complex<double>> values; // point/time/channel
    std::vector<std::string> solutionFamilies; // one per IC
};
// g(Phi^{offset}(x)) on each supplied initial state. Uses registered analytic
// evaluators on their physical/native clock; no integration or Fourier fit.
// Equilibria and nonperiodic trajectories are allowed if their evaluator exists.
AnalyticHistories sampleAnalyticHistories(const std::vector<AnalyticSolutionFamily>&,
    const std::vector<double>&initialStates,std::size_t dimension,
    const std::vector<double>&timeOffsets,const std::vector<Observable>&,
    double initialTime=0,const std::string&family="auto");
}
