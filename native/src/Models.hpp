#pragma once

#include <cstddef>
#include <string>

namespace dynsyslab::core::detail {

void evaluateRHS(const std::string& model, const double* parameters, double time,
                 const double* state, double* output);
void validateState(const std::string& model, const double* state);
std::size_t parameterCount(const std::string& model);
std::size_t stateDimension(const std::string& model);

}  // namespace dynsyslab::core::detail
