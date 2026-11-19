#pragma once
// Classical frequency-domain analysis for the `linear_systems` method: per-channel
// transfer functions and Bode data, and loop robustness (classical and disk margins,
// sensitivity peaks) of the LQR, Kalman filter and LQG loops, with an optional
// loop-transfer-recovery sweep. Appends arrays and metadata to the method's archive.
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/StateSpace.hpp>

#include <string>
#include <vector>

namespace dynsyslab::workflow::postprocess_detail {

void appendTransferFunctions(core::NpzRecord& out, core::json::Value::Object& details,
                             const core::StateSpaceModel& model, const std::vector<std::string>& outputNames,
                             const std::vector<double>& frequencies);

// K-ROM readouts (e.g. state = (r, r-dot)) seen through the measured output: the full
// observability gramian's singular values and directions, each direction's readout
// signature C_r v_k, and (with LQE weights) the Kalman filter's own predicted readout
// standard deviations sqrt(diag(C_r P C_r^T)).
void appendReadoutObservability(core::NpzRecord& out, core::json::Value::Object& details,
                                const core::StateSpaceModel& model, const core::StateSpaceModel& readouts,
                                const std::vector<std::string>& readoutNames, const core::json::Value& settings);

void appendLoopRobustness(core::NpzRecord& out, core::json::Value::Object& details,
                          const core::StateSpaceModel& model, const core::json::Value& settings);

}  // namespace dynsyslab::workflow::postprocess_detail
