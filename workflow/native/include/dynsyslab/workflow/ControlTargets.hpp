#pragma once

#include <dynsyslab/core/LQR.hpp>
#include <dynsyslab/workflow/Json.hpp>

#include <cstddef>

namespace dynsyslab::workflow {

// Resolve a versioned solve-job target into a model-coordinate state function.
// Supported targets are fixed points, kinematic waveforms, archived
// trajectories, and retained canonical targets. Any of them may carry an
// optional "ramp" block (smooth onset from a resting state); "ramp" absent,
// null or false leaves the target unconstrained.
[[nodiscard]] dynsyslab::core::StateTargetFunction resolveControlTarget(
    const json::Value& specification, std::size_t stateDimension);

}  // namespace dynsyslab::workflow
