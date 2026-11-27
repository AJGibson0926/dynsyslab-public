#pragma once
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/HistoryInputRegression.hpp>
#include <dynsyslab/core/HistoryKoopman.hpp>
#include <dynsyslab/core/Integrator.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <filesystem>
#include <functional>

namespace echo_experiment {
namespace core = dynsyslab::core;
using Complex = std::complex<double>;
using Json = core::json::Value;
struct Frozen {
    core::HankelLeftBasis basis;
    core::EigenfunctionFitResult polynomial;
    double sensorStep{}, controlStep{0.01};
    std::size_t stride{20};
    Complex multiplier;
};
struct Measurements { std::vector<double> time, echo, input; };
struct Lifted {
    std::vector<double> time;
    std::vector<std::array<double,2>> xi;
    std::vector<Complex> psi;
    std::vector<core::HistoryInputSample> pairs;
};
Frozen frozen(const core::NpzRecord&, const core::EigenfunctionFitResult&);
Measurements readMeasurements(const std::filesystem::path&);
Lifted lift(const Measurements&, const Frozen&);
using EchoController = std::function<double(double,double)>;
void simulate(const core::SystemConfiguration&, const std::filesystem::path&,
    const std::string&, const Frozen&, double duration, double initialRadius,
    EchoController, Json metadata);
Json controlTrial(const core::SystemConfiguration&, const std::filesystem::path&,
    const std::string&, const Frozen&, const core::HistoryInputMap&,
    const Measurements&, const Lifted&, double feedbackMultiplier);
Json metrics(const Measurements&, const Measurements&, double start);
void saveJson(const std::filesystem::path&,const Json&);
void refineExperiment(const core::SystemConfiguration&,const std::filesystem::path&,
    const Frozen&,const core::NpzRecord&);
void eigenfunctionSweep(const core::SystemConfiguration&,const std::filesystem::path&,
    const Frozen&,const core::NpzRecord&,const std::filesystem::path& inputs = {});
void shortHistoryStudy(const core::SystemConfiguration&,const std::filesystem::path&,
    const Frozen&);
} // namespace echo_experiment
