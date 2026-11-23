#pragma once

#include "CliArguments.hpp"

#include <dynsyslab/workflow/Json.hpp>

namespace dynsyslab::cli {

// Copies the job-group command-line options into the job's "runtime" block.
void applyJobGroupRuntime(dynsyslab::workflow::json::Value& job, const Arguments& args);

// dynsyslab-cli datapack ...: a thin layer over the job module. Execution is executeJob on the
// pack's root job; this command adds pin checks, verification, sealing and export.
int dataPackCommand(const Arguments& args);

}  // namespace dynsyslab::cli
