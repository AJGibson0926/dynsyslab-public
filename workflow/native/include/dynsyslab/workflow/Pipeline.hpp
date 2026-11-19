#pragma once

#include <dynsyslab/core/Execution.hpp>
#include <dynsyslab/workflow/Json.hpp>
#include <dynsyslab/workflow/Run.hpp>

#include <filesystem>
#include <string>
#include <utility>
#include <vector>

namespace dynsyslab::workflow {

// A pipeline job ("kind":"pipeline") is an ordered list of ordinary compute
// jobs ("stages") whose inputs may name earlier stages' published outputs:
//
//   ${id.archive}         the stage's single published archive
//   ${id.archives[N]}     its N-th published archive (zero based)
//   ${id.output:METHOD}   the archive a postprocess stage produced for METHOD
//   ${id.manifest}        its hierarchy/campaign/pipeline manifest
//
// "extends":"base.json" inherits a base pipeline (paths resolve against the
// process working directory, like every other job path). "stages" entries in
// the derived job replace base stages with the same id or are appended, and
// "overrides":{"id":{...}} applies an RFC 7386 JSON merge patch to that stage's
// job. Variants are therefore small, auditable differences from one base.
//
// "parameters":{...} declares named values (merged along the extends chain);
// "${params.NAME}" as a whole string is replaced by that JSON value, so a
// variant can change e.g. SINDy coordinates without restating method arrays.
// A derived stage with "before":"id" is inserted ahead of that base stage.

// Resolve "extends" and "overrides" into one self-contained pipeline job.
[[nodiscard]] json::Value expandPipeline(const json::Value& job);

// Structural validation plus a per-stage plan; stages whose inputs depend on
// earlier outputs are planned when they run.
[[nodiscard]] json::Value planPipeline(const json::Value& job);

[[nodiscard]] RunResult executePipeline(
    const json::Value& job, const std::filesystem::path& archiveRoot,
    const std::filesystem::path& progressPath = {},
    const dynsyslab::core::ExecutionContext* execution = nullptr);

// Command-line style inputs, so no file needs editing to vary a run:
//   params: NAME=VALUE  -> pipeline "parameters" (merged before expansion)
//   sets:   PATH=VALUE  -> any field by dotted path; numeric segments index
//                          arrays. For a pipeline the first segment is a stage
//                          id and the rest addresses that stage's job, e.g.
//                          "control.controller.settings.Q=1000".
// VALUE is parsed as JSON, falling back to a plain string. A pipeline is
// returned fully expanded (no "extends"), with the inputs recorded under
// "command_inputs" for the manifest.
[[nodiscard]] json::Value applyJobInputs(
    json::Value job, const std::vector<std::pair<std::string, std::string>>& params,
    const std::vector<std::pair<std::string, std::string>>& sets);

// RFC 7386 merge patch (objects merge recursively, null deletes, anything else
// replaces). Exposed for tests and for callers composing job variants.
[[nodiscard]] json::Value mergePatch(json::Value target, const json::Value& patch);

}  // namespace dynsyslab::workflow
