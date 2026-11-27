#pragma once

// Numerical settings of an ODE `solve` job (docs/architecture/NATIVE_SOLVER.md):
// method, output clock, guards and step budget, failure policy, post-solve resampling,
// checkpoints and pause.

#include <dynsyslab/core/Integrator.hpp>
#include <dynsyslab/core/TimeSeries.hpp>
#include <dynsyslab/workflow/Json.hpp>

#include <filesystem>
#include <string>
#include <vector>

namespace dynsyslab::workflow::detail {

enum class FailurePolicy { error, truncate, reject };

struct ResampleSettings {
    bool enabled{false};
    core::Interpolation method{core::Interpolation::pchip};
    std::vector<double> times;   // explicit, strictly increasing, inside the solved span
};

struct CheckpointSettings {
    bool enabled{false};
    double interval{0.0};                // model time between partial archives
    std::filesystem::path pauseFile;     // pause when this file exists at a checkpoint
    double wallInterval{0.0};            // wall-clock seconds between partial archives
    double deadlineUnix{0.0};            // checkpoint and pause at this UTC time (s since epoch)
};

struct SolveSettings {
    core::IntegrationOptions integration;
    FailurePolicy failure{FailurePolicy::error};
    ResampleSettings resample;
    CheckpointSettings checkpoints;
};

// Reads job.time, job.solver, job.interpolation and job.history. Time values are in the
// integration coordinates (the system's model clock, or seconds when
// time.coordinate_system is "dimensional").
[[nodiscard]] SolveSettings solveSettings(const json::Value& job);

[[nodiscard]] const char* failurePolicyName(FailurePolicy policy);

}  // namespace dynsyslab::workflow::detail
