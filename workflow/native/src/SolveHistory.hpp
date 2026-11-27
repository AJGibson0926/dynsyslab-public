#pragma once

// History options of an ODE `solve` job beyond checkpoints (docs/architecture/NATIVE_SOLVER.md):
//   history.controller_call_log (bool), history.controller_call_log_limit (entries, 0 = all):
//     every controller evaluation is archived per trajectory as controller_call_time,
//     controller_call_state, controller_call_input, controller_call_kind (0 right-hand side,
//     1 sampled-control tick, 2 output sample) and controller_call_status (1 accepted step,
//     0 rejected trial step, -1 unknown).
//   history.keep_partial_history_on_failure (bool): with solver.on_failure "error", a failed
//     solve still fails, but first writes each started trajectory's computed samples to
//     <archive>.partial.npz (schema dynsyslab.solution-checkpoint/v1, kept_on_failure true).

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/ControllerCallLog.hpp>
#include <dynsyslab/core/Integrator.hpp>
#include <dynsyslab/workflow/Json.hpp>

#include <cstddef>
#include <filesystem>
#include <functional>
#include <string>

namespace dynsyslab::workflow::detail {

// Reads job.history and attaches `log` and the keep-partial flag to `options`.
void applyHistoryOptions(const json::Value& job, core::IntegrationOptions& options,
                         core::ControllerCallLog& log, std::size_t stateDimension,
                         std::size_t inputDimension);

// Archives one trajectory's controller calls (no-op when the log is not attached).
void appendControllerCalls(core::NpzRecord& archive, json::Value::Object& metadata,
                           const core::IntegrationOptions& options, std::size_t trajectory);

// Runs the integration; on core::IntegrationFailure writes the partial archives and rethrows
// the failure as std::runtime_error naming them.
[[nodiscard]] core::NativeSolution integrateKeepingHistory(
    const std::string& label, const std::function<std::filesystem::path(std::size_t)>& archivePath,
    const std::function<core::NativeSolution()>& run);

}  // namespace dynsyslab::workflow::detail
