#pragma once
// Experiment stage kinds beyond plain jobs (EXPERIMENTS.md, "Sweeps and selection"):
//
//   sweep   {"id": "delays", "sweep": {"variables": {"delays": [50, 100], "rank": [2, 4]},
//            "mode": "product" | "zip", "stages": [{"id": "basis", "job": {...}}, ...]}}
//           copies its block of stages once per point (product: first variable slowest), as
//           DELAYS_K_BASIS, with ${sweep.NAME} replaced by the point's value (a whole string
//           becomes the JSON value) and references to the block's own stages renamed.
//   select  {"id": "best", "select": {"sweep": "delays", "stage": "score", "object": "report",
//            "metric": "validation.crossings", "goal": "minimize", "reduce": "mean"}}
//           scores each point's STAGE output by a metadata path (or "array:NAME") of a
//           validation section and stands for the winner; ${best.archive} is its archive (or, with
//           "stands_for": INNER, the archive of that other stage of the winning point) and
//           ${best.manifest} the selection record (every candidate, its score, the rejects).
#include <dynsyslab/workflow/Json.hpp>
#include <dynsyslab/workflow/Run.hpp>

#include <filesystem>

namespace dynsyslab::workflow::experiment_detail {

// Replaces sweep stages by their copies and select stages by "select" jobs over the candidates;
// records the points under "sweeps".
void expandSweepsAndSelects(json::Value& expanded);

// Runs a resolved "select" job and writes its selection record to `record`.
[[nodiscard]] RunResult runSelect(const json::Value& job, const std::filesystem::path& record);

}  // namespace dynsyslab::workflow::experiment_detail
