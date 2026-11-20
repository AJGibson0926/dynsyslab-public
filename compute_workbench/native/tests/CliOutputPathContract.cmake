# Output paths given on the command line: `run --progress progress.jsonl` (bare name, resolved
# against the working directory), a relative sub-folder, an absolute path, and an uncreatable
# parent (a JSON error naming the option, before any work). Also DynSysLabViewerCli --output plot.png.
file(REMOVE_RECURSE "${ROOT}")
file(MAKE_DIRECTORY "${ROOT}/cwd")
set(archives "${ROOT}/archives")

function(run_job progress)
  execute_process(COMMAND "${CLI}" run "${JOB}" --archive-root "${archives}" --progress "${progress}"
    WORKING_DIRECTORY "${ROOT}/cwd" RESULT_VARIABLE code OUTPUT_VARIABLE out ERROR_VARIABLE err)
  set(code "${code}" PARENT_SCOPE)
  set(out "${out}" PARENT_SCOPE)
  set(err "${err}" PARENT_SCOPE)
endfunction()

run_job("progress.jsonl")
if(NOT code EQUAL 0 OR NOT EXISTS "${ROOT}/cwd/progress.jsonl")
  message(FATAL_ERROR "bare --progress file name failed (${code}): ${err}")
endif()
run_job("sub/deeper/progress.jsonl")
if(NOT code EQUAL 0 OR NOT EXISTS "${ROOT}/cwd/sub/deeper/progress.jsonl")
  message(FATAL_ERROR "relative --progress sub-folder failed (${code}): ${err}")
endif()
run_job("${ROOT}/absolute/progress.jsonl")
if(NOT code EQUAL 0 OR NOT EXISTS "${ROOT}/absolute/progress.jsonl")
  message(FATAL_ERROR "absolute --progress failed (${code}): ${err}")
endif()

file(GLOB before "${archives}/solutions/*/*.npz")
list(LENGTH before archivesBefore)
file(WRITE "${ROOT}/cwd/blocker" "a file where a folder is needed")
run_job("blocker/progress.jsonl")
if(code EQUAL 0)
  message(FATAL_ERROR "an uncreatable --progress parent was accepted")
endif()
string(JSON schema GET "${err}" schema)
string(JSON option GET "${err}" option)
string(JSON path GET "${err}" path)
if(NOT schema STREQUAL "dynsyslab.cli-error/v1" OR NOT option STREQUAL "--progress" OR NOT path MATCHES "blocker")
  message(FATAL_ERROR "path error is not a dynsyslab.cli-error/v1 naming --progress and the path: ${err}")
endif()
file(GLOB after "${archives}/solutions/*/*.npz")
list(LENGTH after archivesAfter)
if(NOT archivesAfter EQUAL archivesBefore)
  message(FATAL_ERROR "a bad --progress path still ran the job")
endif()

list(GET after 0 archive)
execute_process(COMMAND "${VIEWER}" time-series --source "${archive}" --curve "values=state"
    --output plot.png --width 600 --height 400
  WORKING_DIRECTORY "${ROOT}/cwd" RESULT_VARIABLE code OUTPUT_VARIABLE out ERROR_VARIABLE err)
if(NOT code EQUAL 0 OR NOT EXISTS "${ROOT}/cwd/plot.png")
  message(FATAL_ERROR "DynSysLabViewerCli --output plot.png failed (${code}): ${err}")
endif()
execute_process(COMMAND "${VIEWER}" time-series --source "${archive}" --curve "values=state"
    --output "blocker/plot.png"
  WORKING_DIRECTORY "${ROOT}/cwd" RESULT_VARIABLE code OUTPUT_VARIABLE out ERROR_VARIABLE err)
string(JSON option GET "${err}" option)
if(code EQUAL 0 OR NOT option STREQUAL "--output")
  message(FATAL_ERROR "Viewer path error does not name --output: ${err}")
endif()
