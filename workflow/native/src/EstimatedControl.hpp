#pragma once
#include <dynsyslab/workflow/ControlFactory.hpp>
#include <functional>
namespace dynsyslab::workflow::detail {
ControllerSetup buildEstimatedController(const json::Value& controller,
    const json::Value& solver, const core::ODESystem& system,
    std::size_t trajectoryCount, double defaultSampleTime,
    std::function<double(std::size_t,double,const double*)> appliedOverride = {});
}
