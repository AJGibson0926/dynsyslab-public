#pragma once

#include <dynsyslab/core/Execution.hpp>
#include <dynsyslab/workflow/Json.hpp>
#include <dynsyslab/workflow/Run.hpp>

#include <cstddef>
#include <filesystem>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace dynsyslab::workflow {

// A job group ("kind":"job_group") is an ordinary compute job whose work is other job files:
//
//   {"schema":"dynsyslab.compute-job/v1","kind":"job_group","label":"study",
//    "parameters":{"Q":5000},
//    "cleanup":{"enabled":true,"min_free_gb":6},
//    "jobs":[{"id":"train","job":"10-data/train.job.json"},
//            {"id":"krom","job":"20-krom/krom.job.json",
//             "inputs":{"trajectory":"${train.outputs.trajectory}"},
//             "parameters":{"degree":4}}],
//    "outputs":{"model":"${krom.outputs.model}"}}
//
// Child paths are relative to the file that names them. A child may be another group. Inside a
// child, ${inputs.NAME} and ${params.NAME} are replaced; a child declares its archives with
// "outputs":{"NAME":"archive"|"archives[N]"|"output:METHOD"|"metadata:KEY"}. Children run in dependency order
// (declared order among independent children); a failed child blocks only its dependents.
// The tree is flattened to leaf jobs, each of which runs through executeJob.

struct JobGroupLeaf {
    std::string id;            // path of group entry ids, e.g. "10-data/train"
    std::filesystem::path file;
    json::Value job;           // inputs bound to @{leaf|selector} tokens
    std::set<std::string> dependencies;
    bool keep{};               // survives cleanup
    double estimatedBytes{};   // "resources":{"peak_bytes":N}, for the disk preflight
};

struct JobGroupGraph {
    std::vector<JobGroupLeaf> leaves;
    std::map<std::string, std::string> outputs;  // root outputs: name -> token
};

// Records the file a job was read from, so a group can resolve children relative to it.
// Any other kind of job is returned unchanged.
[[nodiscard]] json::Value bindJobFile(json::Value job, const std::filesystem::path& file);

[[nodiscard]] JobGroupGraph loadJobGroup(const json::Value& job);
[[nodiscard]] json::Value planJobGroup(const json::Value& job);
[[nodiscard]] RunResult executeJobGroup(
    const json::Value& job, const std::filesystem::path& archiveRoot,
    const std::filesystem::path& progressPath = {},
    const dynsyslab::core::ExecutionContext* execution = nullptr);

// The folder holding a group's ledger, leaf archives and figures.
[[nodiscard]] std::filesystem::path jobGroupRoot(const json::Value& job,
                                                 const std::filesystem::path& archiveRoot);

}  // namespace dynsyslab::workflow
