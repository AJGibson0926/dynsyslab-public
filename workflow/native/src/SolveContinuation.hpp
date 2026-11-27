#pragma once

// Continuing a solution: extend an archive past its stop time, or resume a paused solve.
// Archives are immutable, so a continuation is a new child archive whose metadata links its
// parent (path and SHA-256) and the parent time it starts from.

#include <dynsyslab/core/Integrator.hpp>
#include <dynsyslab/workflow/Json.hpp>

#include <filesystem>
#include <string>
#include <vector>

namespace dynsyslab::workflow::detail {

struct ContinuationPoint {
    double time{};
    std::vector<double> state;        // last computed (finite) row of the parent
    std::size_t completedSamples{};
};

[[nodiscard]] ContinuationPoint continuationPoint(const std::filesystem::path& archive);

// initial_conditions {mode: "continue_from", archive | archives}.
[[nodiscard]] std::vector<std::string> continuationArchives(const json::Value& initialConditions);

// The `continuation` metadata block of trajectory `index` (empty object when the job does
// not continue an archive): parent_archive, parent_archive_sha256, parent_time.
[[nodiscard]] json::Value continuationMetadata(const json::Value& job, std::size_t index);

// Throws unless every parent ends at time.start.
void requireContinuationStart(const json::Value& job, double start);

// A complete solve job that continues `parent` from its last computed sample to `stop`,
// built from the parent's archived system, controller, solver and time configuration.
[[nodiscard]] json::Value extensionJob(const std::filesystem::path& parent, double stop);

// When trajectories were rejected by guards or a pause stopped the solve, writes
// <manifest>.outcome.json (dynsyslab.solve-outcome/v1): the accepted-member mask, each
// rejected trajectory's reason and, after a pause, `resume_jobs` -- a continuation of the
// paused archive and a solve of the trajectories that had not started. Running those jobs
// resumes the solve; the parent archives stay unchanged.
void writeSolveOutcome(const std::filesystem::path& manifest, const json::Value& job,
                       const core::NativeSolution& solution,
                       const std::vector<std::filesystem::path>& archives,
                       const std::vector<std::size_t>& kept, const std::vector<double>& initials);

}  // namespace dynsyslab::workflow::detail
