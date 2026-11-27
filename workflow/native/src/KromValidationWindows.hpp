#pragma once
// koopman_validation windows, binned defects and psi history (KromValidationWindows.cpp).
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/KoopmanModel.hpp>
#include <dynsyslab/workflow/Json.hpp>

#include <vector>

namespace dynsyslab::workflow::postprocess_detail {

void appendValidationWindows(dynsyslab::core::NpzRecord& archive, json::Value& meta,
                             const dynsyslab::core::KoopmanModel& model, const std::vector<double>& time,
                             const std::vector<double>& values, std::size_t samples, std::size_t dimension,
                             const json::Value& settings);

}  // namespace dynsyslab::workflow::postprocess_detail
