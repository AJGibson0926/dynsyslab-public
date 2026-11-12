#pragma once
// Isolated experiment composition; no primary-pipeline changes.
#include "Experiment.hpp"
#include <dynsyslab/core/KoopmanValidation.hpp>
#include <dynsyslab/core/LQR.hpp>
namespace echo_multi {
using namespace echo_experiment;
namespace fs = std::filesystem;
struct Coordinates { std::vector<double> time, x, input; };
Coordinates project(const Measurements&, const Frozen&, double begin=0, double end=1e300);
Json freeErrors(const core::KoopmanModel&, const std::vector<Coordinates>&, double, double);
void runControls(const fs::path& root,const fs::path& output,const std::string& window,
    const Frozen&,const core::KoopmanModel&,const core::SystemConfiguration&,bool joint=false);
void auditResults(const fs::path&,const fs::path&);
using ModeSamples = std::vector<std::vector<core::HistoryInputSample>>;
struct JointInputFit {
    core::BilinearControlMap map;
    std::vector<Complex> coefficients; // output x [1, Re z1, Im z1, ...], discrete
    double condition{};
    std::size_t rank{},columns{},samples{};
};
JointInputFit fitJointInput(const ModeSamples&,const std::vector<std::size_t>&,
    const std::vector<Complex>& eigenvalues,double step=.01,double ridge=1e-8);
Json scoreJointInput(const JointInputFit&,const ModeSamples&,const std::vector<std::size_t>&,
    const std::vector<Complex>& eigenvalues,double step=.01);
Json jointInputTests();
Json runModelTrial(const fs::path&,const std::string&,const Frozen&,
    const core::KoopmanModel&,const core::SystemConfiguration&,const Measurements&,
    const Coordinates&,const std::vector<double>& scales,double penalty);
}
