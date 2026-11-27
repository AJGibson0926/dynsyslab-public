#pragma once
// dynsyslab-cli queue | campaign | benchmark | template (docs/architecture/JOB_QUEUE.md).
//   queue submit <job.json> [--retry-limit N] [--after ID[:CONDITION]]... [--hold]
//                [--reuse identical] [--max-concurrent N]
//   queue list | status ID | cancel ID | pause ID | resume ID
//   queue tray | release ID | remove ID | watchdog [--grace-seconds S] | recover AGE | resources
//   queue run-next [--max-threads N] [--max-memory-mib M]
//   queue serve [--workers N] [--max-running N] [--pool-threads N] [--pool-memory-mib M]
//               [--max-threads N] [--max-memory-mib M] [--poll-seconds S] [--max-jobs N]
//               [--idle-exit S] [--grace-seconds S] [--command-template "CMD {job_file}"]
//   campaign list | status ID | validate ID | repair ID | profiles
//   benchmark [run] [--iterations N] [--warmup N] [--case NAME]... [--output FILE]
//   benchmark gate <reference.json> <candidate.json> [--required-speedup X]
//             [--require CASE=X]... [--any-runtime]
//   template solve|postprocess|estimate [--output FILE]
// Options after the command are read from the positional tokens CliArguments leaves.
#include "CliArguments.hpp"

#include <dynsyslab/core/Execution.hpp>

namespace dynsyslab::cli {

int queueCommand(const Arguments& args);
int campaignCommand(const Arguments& args);
int benchmarkCommand(const Arguments& args);
int templateCommand(const Arguments& args);
// --max-threads N and --max-memory-mib M among the positional tokens (0 = unlimited).
[[nodiscard]] dynsyslab::core::ExecutionLimits executionLimits(const Arguments& args);

}  // namespace dynsyslab::cli
