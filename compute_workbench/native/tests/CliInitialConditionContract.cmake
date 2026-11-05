execute_process(
  COMMAND "${CLI}" generate-ics "${JOB}"
  RESULT_VARIABLE status
  OUTPUT_VARIABLE output
  ERROR_VARIABLE error)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "generate-ics failed: ${error}")
endif()
string(JSON schema GET "${output}" schema)
string(JSON count GET "${output}" trajectory_count)
string(JSON first_x GET "${output}" values 0 0)
if(NOT schema STREQUAL "dynsyslab.initial-conditions/v1")
  message(FATAL_ERROR "unexpected IC response schema: ${schema}")
endif()
if(NOT count EQUAL 1 OR NOT first_x EQUAL -1)
  message(FATAL_ERROR "unexpected generated IC response: ${output}")
endif()
