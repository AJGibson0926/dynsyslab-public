#pragma once

#include <dynsyslab/core/LQR.hpp>
#include <dynsyslab/workflow/Json.hpp>

#include <cstddef>
#include <string>
#include <vector>

namespace dynsyslab::workflow {

// Resolve a versioned solve-job target into a model-coordinate state function.
// Supported targets are fixed points, kinematic waveforms, archived
// trajectories, and retained canonical targets. Any of them may carry an
// optional "ramp" block (smooth onset from a resting state); "ramp" absent,
// null or false leaves the target unconstrained.
// Kinds: constant | fixed_point, waveform, trajectory (archived; "interpolation" linear or
// pchip), expression (one analytic expression of t per component), rpe_canonical.
[[nodiscard]] dynsyslab::core::StateTargetFunction resolveControlTarget(
    const json::Value& specification, std::size_t stateDimension);

// Bounds that contain the target on [start, stop]: analytic for harmonic targets (Python
// conservative_state_envelope), otherwise sampled at 4001 times.
struct TargetEnvelope {
    std::vector<double> lower, upper;
    std::string method;
};
[[nodiscard]] TargetEnvelope controlTargetEnvelope(const json::Value& specification,
                                                   std::size_t stateDimension, double start, double stop);

}  // namespace dynsyslab::workflow
