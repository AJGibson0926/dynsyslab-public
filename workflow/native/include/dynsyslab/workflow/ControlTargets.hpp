#pragma once

#include <dynsyslab/core/LQR.hpp>
#include <dynsyslab/workflow/Json.hpp>

#include <cstddef>

namespace dynsyslab::workflow {

// Resolve a versioned solve-job target into a model-coordinate state function.
// Supported targets are fixed points, kinematic waveforms, archived
// trajectories, and retained canonical targets.
[[nodiscard]] dynsyslab::core::StateTargetFunction resolveControlTarget(
    const json::Value& specification, std::size_t stateDimension);

}  // namespace dynsyslab::workflow
