cmake_policy(SET CMP0054 NEW)

if(NOT EXISTS "${CLI}")
  message(FATAL_ERROR "compiled CLI is missing: ${CLI}")
endif()
if(NOT EXISTS "${JOB}")
  message(FATAL_ERROR "queue test job is missing: ${JOB}")
endif()
file(REMOVE_RECURSE "${ROOT}")

execute_process(COMMAND "${CLI}" queue submit "${JOB}" --archive-root "${ROOT}"
  RESULT_VARIABLE status OUTPUT_VARIABLE submitted ERROR_VARIABLE error)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "queue submit failed: ${error}")
endif()
string(JSON schema GET "${submitted}" schema)
string(JSON job_id GET "${submitted}" payload job_id)
if(NOT schema STREQUAL "dynsyslab.job-queue/v1" OR job_id STREQUAL "")
  message(FATAL_ERROR "queue submit JSON contract failed")
endif()

execute_process(COMMAND "${CLI}" queue status "${job_id}" --archive-root "${ROOT}"
  RESULT_VARIABLE status OUTPUT_VARIABLE queued ERROR_VARIABLE error)
string(JSON job_status GET "${queued}" payload status)
if(NOT status EQUAL 0 OR NOT job_status STREQUAL "queued")
  message(FATAL_ERROR "queue status JSON contract failed: ${error}")
endif()

execute_process(COMMAND "${CLI}" queue cancel "${job_id}" --archive-root "${ROOT}"
  RESULT_VARIABLE status OUTPUT_VARIABLE cancelled ERROR_VARIABLE error)
string(JSON job_status GET "${cancelled}" payload status)
if(NOT status EQUAL 0 OR NOT job_status STREQUAL "cancelled")
  message(FATAL_ERROR "queue cancel JSON contract failed: ${error}")
endif()

execute_process(COMMAND "${CLI}" queue submit "${JOB}" --archive-root "${ROOT}"
  RESULT_VARIABLE status OUTPUT_VARIABLE submitted ERROR_VARIABLE error)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "second queue submit failed: ${error}")
endif()
execute_process(COMMAND "${CLI}" queue run-next --archive-root "${ROOT}"
  RESULT_VARIABLE status OUTPUT_VARIABLE completed ERROR_VARIABLE error)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "queue run-next failed: ${error}")
endif()
string(JSON queue_schema GET "${completed}" schema)
string(JSON result_schema GET "${completed}" payload result schema)
string(JSON result_kind GET "${completed}" payload result kind)
string(JSON archive_count LENGTH "${completed}" payload result archives)
if(NOT queue_schema STREQUAL "dynsyslab.job-queue/v1" OR
   NOT result_schema STREQUAL "dynsyslab.workflow-result/v1" OR
   NOT result_kind STREQUAL "solve" OR archive_count LESS 1)
  message(FATAL_ERROR "queue execution result JSON contract failed")
endif()

# Phase C: retry flag, side tray, templates, campaign profiles, benchmark cases and gates.
execute_process(COMMAND "${CLI}" queue submit "${JOB}" --hold --retry-limit 2 --archive-root "${ROOT}"
  RESULT_VARIABLE status OUTPUT_VARIABLE held ERROR_VARIABLE error)
string(JSON held_id GET "${held}" payload job_id)
execute_process(COMMAND "${CLI}" queue status "${held_id}" --archive-root "${ROOT}"
  RESULT_VARIABLE status OUTPUT_VARIABLE reserved ERROR_VARIABLE error)
string(JSON held_status GET "${reserved}" payload status)
string(JSON held_retry GET "${reserved}" payload retry_limit)
if(NOT status EQUAL 0 OR NOT held_status STREQUAL "reserved" OR NOT held_retry EQUAL 2)
  message(FATAL_ERROR "queue side tray or retry flag failed: ${reserved} ${error}")
endif()
execute_process(COMMAND "${CLI}" queue remove "${held_id}" --archive-root "${ROOT}"
  RESULT_VARIABLE status ERROR_VARIABLE error)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "queue remove failed: ${error}")
endif()
execute_process(COMMAND "${CLI}" template solve
  RESULT_VARIABLE status OUTPUT_VARIABLE template ERROR_VARIABLE error)
string(JSON template_kind GET "${template}" kind)
if(NOT status EQUAL 0 OR NOT template_kind STREQUAL "solve")
  message(FATAL_ERROR "template solve failed: ${error}")
endif()
execute_process(COMMAND "${CLI}" campaign profiles
  RESULT_VARIABLE status OUTPUT_VARIABLE profiles ERROR_VARIABLE error)
string(JSON low_points GET "${profiles}" quality low points)
if(NOT status EQUAL 0 OR NOT low_points EQUAL 31)
  message(FATAL_ERROR "campaign profiles failed: ${error}")
endif()
execute_process(COMMAND "${CLI}" benchmark --case ode_rk4_duffing --iterations 1 --warmup 0
  --output "${ROOT}/benchmark.json" RESULT_VARIABLE status ERROR_VARIABLE error)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "benchmark failed: ${error}")
endif()
execute_process(COMMAND "${CLI}" benchmark gate "${ROOT}/benchmark.json" "${ROOT}/benchmark.json"
  --required-speedup 0.5 RESULT_VARIABLE status OUTPUT_VARIABLE gate ERROR_VARIABLE error)
string(JSON gate_passed GET "${gate}" passed)
if(NOT status EQUAL 0 OR NOT gate_passed)
  message(FATAL_ERROR "benchmark gate failed: ${gate} ${error}")
endif()
file(REMOVE_RECURSE "${ROOT}")
