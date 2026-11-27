#pragma once
// Non-executing K-ROM training-domain expansion plans (`dynsyslab-cli krom-training-plan`).
//
// Request (dynsyslab.krom-training-plan-request/v1):
//   "krom": PATH (its training domain is the current domain) | "current_domain": {minimum, maximum}
//   "system": {name, preset, parameters}       default: the K-ROM's system
//   "requirements": [{"name", "domain": {minimum, maximum}},
//                    {"name", "target": TARGET, "probes": 33, "start": 0, "stop": 2 pi}, ...]
//       a target requirement covers the target's envelope (analytic for harmonic targets) and adds
//       `probes` target states, equally spaced in t on [start, stop), that must survive screening
//   "preset": "rpe_targets" (with optional "targets": [names]) adds every canonical RPE target
//   "points_per_axis": [...], "margin_fraction": 0.1, "anchors": [[...]],
//   "positive_axes": [i, ...]  refuse a plan whose planned minimum is <= 0 on these axes (r > 0)
//   "time": {"stop": 25, "output_step": 0.001}, "solver": {...}, "screening": {...}
// Result (dynsyslab.krom-training-plan/v1): current, required (the union of every requirement),
// planned domains; anchors and named probes; initial-condition count; plan_fingerprint_sha256;
// and "solve_job", a solve job whose explicit initial conditions are the planned grid plus
// anchors and probes, with the krom_training_domain preflight that screens them all.
#include <dynsyslab/workflow/Json.hpp>

namespace dynsyslab::workflow {

[[nodiscard]] json::Value designKromTrainingPlan(const json::Value& request);

}  // namespace dynsyslab::workflow
