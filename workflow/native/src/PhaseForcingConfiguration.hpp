#pragma once
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/System.hpp>
namespace dynsyslab::workflow::detail {
dynsyslab::core::SinusoidalDriver
phaseDriverFromController(const dynsyslab::core::json::Value &controller,
                          const std::string &evaluation, const dynsyslab::core::ODESystem &system);
}
