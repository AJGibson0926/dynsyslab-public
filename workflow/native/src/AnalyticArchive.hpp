#pragma once
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/System.hpp>
#include <dynsyslab/core/TimeSeries.hpp>
namespace dynsyslab::workflow::detail {
void appendAnalyticComparison(dynsyslab::core::NpzRecord &archive,
                              const dynsyslab::core::SystemDescriptor &descriptor,
                              const dynsyslab::core::SystemConfiguration &configuration,
                              const std::vector<dynsyslab::core::AnalyticSolutionFamily> &families,
                              const std::vector<double> &time, const std::vector<double> &values,
                              const std::vector<double> &initial, double initialTime,
                              const dynsyslab::core::json::Value &job);
void appendAnalyticComparison(dynsyslab::core::NpzRecord &archive,
                              const dynsyslab::core::ODESystem &system,
                              const dynsyslab::core::TimeSeries &series,
                              const std::vector<double> &initial, double initialTime,
                              const dynsyslab::core::json::Value &job);
} // namespace dynsyslab::workflow::detail
