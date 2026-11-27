#pragma once
// control_score setting "estimation": {"components":[0,1], "innovation_lags":[0.01,0.1,0.5],
// "shadow_force":true}. Per window, on the estimator clock:
//   estimate_error_rms_K, estimate_bias_K   state estimate minus the solved state, component K
//   state_std_K, estimate_relative_error_K  the state's standard deviation and rms error / std
//   innovation_relative                     innovation rms / measurement standard deviation
//   innovation_autocorrelation_LAG          innovation whiteness at a time lag (model time units)
//   shadow_force_*                          (shadow_force) the archived commanded_input of a shadow
//                                           controller against the applied_input: error rms, rms of
//                                           the applied drive, relative error, correlation
#include "PostprocessInternal.hpp"

#include <dynsyslab/core/Archive.hpp>

#include <map>
#include <string>
#include <vector>

namespace dynsyslab::workflow::postprocess_detail {

[[nodiscard]] std::vector<std::string> estimationMetricNames(const json::Value& estimation);
[[nodiscard]] std::map<std::string, double> scoreEstimation(const core::NpzRecord& record, double start, double stop,
                                                            const json::Value& estimation);

}  // namespace dynsyslab::workflow::postprocess_detail
