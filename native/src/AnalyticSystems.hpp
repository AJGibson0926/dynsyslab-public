#pragma once
#include <dynsyslab/core/System.hpp>
namespace dynsyslab::core::detail {
std::vector<AnalyticSolutionFamily> integrableAnalyticSolutions(const SystemConfiguration &);
std::vector<AnalyticSolutionFamily> todaAnalyticSolutions(const SystemConfiguration &);
std::vector<AnalyticSolutionFamily> triangularAnalyticSolutions(const SystemConfiguration &);
std::vector<AnalyticSolutionFamily> henonHeilesAnalyticSolutions(const SystemConfiguration &);
std::vector<AnalyticSolutionFamily> kuramotoPairAnalyticSolutions(const SystemConfiguration &);
std::vector<AnalyticSolutionFamily> neuralLimitAnalyticSolutions(const SystemConfiguration &);
std::vector<AnalyticSolutionFamily> springLimitAnalyticSolutions(const SystemConfiguration &);
std::vector<AnalyticSolutionFamily> bubbleLimitAnalyticSolutions(const SystemConfiguration &);
std::vector<AnalyticSolutionFamily> oscillatorAnalyticSolutions(const SystemConfiguration &);
std::vector<AnalyticSolutionFamily> doublePendulumAnalyticSolutions(const SystemConfiguration &);
std::vector<AnalyticSolutionFamily> rpeAnalyticSolutions(const SystemConfiguration &);
std::vector<AnalyticSolutionFamily> otherAnalyticSolutions(const SystemDescriptor &,
                                                           const SystemConfiguration &);
std::function<AnalyticState(double)>
linearOscillator(double damping, double stiffness, const AnalyticState &initial, double initialTime,
                 SinusoidalDriver driver = {}, double gain = 1);
AnalyticTrajectory pendulumTrajectory(double frequencySquared, const AnalyticState &initial,
                                      double initialTime);
} // namespace dynsyslab::core::detail
