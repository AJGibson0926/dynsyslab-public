#pragma once
// Pieces of buildController (ControlFactory.cpp) shared with the Koopman controllers
// (ControlFactoryKoopman.cpp).
#include <dynsyslab/workflow/ControlFactory.hpp>
#include <dynsyslab/core/ClassicalControl.hpp>
#include <dynsyslab/core/ODE.hpp>

#include <cstddef>
#include <string>
#include <vector>

namespace dynsyslab::workflow::control_factory_detail {

std::vector<double> numbers(const json::Value& value);
// A penalty given as a scalar, a diagonal or a full dimension x dimension matrix.
std::vector<double> squareMatrix(const json::Value& source, std::size_t dimension, const char* name);
std::size_t positiveCount(const json::Value& source, std::size_t fallback, const char* name);

struct Clock {
    dynsyslab::core::ControlEvaluationMode mode{
        dynsyslab::core::ControlEvaluationMode::continuousRhs};
    std::string label{"continuous_rhs"};
    double sampleTime{};
};
// The controller's evaluation clock from the solver's control_evaluation and sample times.
Clock resolveClock(const std::string& kind, const json::Value& solver,
                   const json::Value& settings, double fallback, bool sampledRequired);

// Multi-input lqr/mpc/pid with an explicit n x m B and additive actuation; false when the
// controller has one input (ControlFactoryClassical.cpp).
bool buildMultiInputController(const std::string& kind, const json::Value& settings, const json::Value& solver,
                               const dynsyslab::core::ODESystem& system, std::size_t trajectoryCount,
                               double defaultSampleTime, const dynsyslab::core::StateTargetFunction& target,
                               ControllerSetup& result);
// A user state-dependent-coefficient factor: rows x columns expressions of the state names.
dynsyslab::core::StateMatrixFunction sdcFactor(const json::Value& entries, std::size_t rows, std::size_t columns,
                                               const dynsyslab::core::ODESystem& system, const char* name);

// KLQR and KMPC on a K-ROM (settings.krom): loads it, binds its history transform and input maps,
// resolves the target in the K-ROM's coordinates, and sets result's policy and evidence.
void buildKoopmanControl(const std::string& kind, const json::Value& settings, const json::Value& targetSpec,
                         const json::Value& solver, const dynsyslab::core::ODESystem& system,
                         std::size_t trajectoryCount, double defaultSampleTime, std::size_t states,
                         bool observableTarget, dynsyslab::core::StateTargetFunction& target,
                         bool& historyClock, ControllerSetup& result);

}  // namespace dynsyslab::workflow::control_factory_detail
