# Exercise the reusable observation stages through the installed native CLI.
# Every artifact lives in the isolated package-smoke directory.
set(_history_job "${STAGE_DIR}/history-smoke-job.json")
file(WRITE "${_history_job}" [=[{
  "schema":"dynsyslab.compute-job/v1","kind":"solve","label":"history-package-source",
  "system":{"name":"duffing","preset":"double_well_reference","parameters":{}},
  "time":{"start":0,"stop":0.4,"output_step":0.01,"coordinate_system":"automatic"},
  "initial_conditions":{"mode":"explicit","count":1,"values":[[0.5,0.1]]},
  "controller":{"kind":"none","settings":{}},
  "solver":{"method":"RK45","rtol":1e-8,"atol":1e-10},
  "archive":{"enabled":true,"data_role":"free_training"}
}]=])
execute_process(COMMAND "${_bin}/dynsyslab-cli.exe" run "${_history_job}"
  --archive-root "${STAGE_DIR}/history-smoke-archives"
  RESULT_VARIABLE _history_status OUTPUT_VARIABLE _history_result ERROR_VARIABLE _history_error)
if(NOT _history_status EQUAL 0)
  message(FATAL_ERROR "packaged history source failed: ${_history_error}")
endif()
require_run_result("${_history_result}" "solve" "packaged history source")
string(JSON _history_source GET "${_history_result}" archives 0)
file(TO_CMAKE_PATH "${_history_source}" _history_source)
foreach(_history_method history_basis history_projection)
  set(_history_settings "\"observable\":\"x\",\"delays\":8,\"hankel_rank\":2,\"coordinate_rank\":2")
  if(_history_method STREQUAL "history_projection")
    string(APPEND _history_settings ",\"transform\":\"${_history_basis}\"")
  endif()
  file(WRITE "${_history_job}" "{\"schema\":\"dynsyslab.compute-job/v1\",\"kind\":\"postprocess\","
    "\"label\":\"history-package-smoke\",\"sources\":[\"${_history_source}\"],"
    "\"methods\":[{\"name\":\"${_history_method}\",\"settings\":{${_history_settings}}}],"
    "\"archive\":{\"enabled\":true}}")
  execute_process(COMMAND "${_bin}/dynsyslab-cli.exe" run "${_history_job}"
    --archive-root "${STAGE_DIR}/history-smoke-archives"
    RESULT_VARIABLE _history_status OUTPUT_VARIABLE _history_result ERROR_VARIABLE _history_error)
  if(NOT _history_status EQUAL 0)
    message(FATAL_ERROR "packaged ${_history_method} failed: ${_history_error}")
  endif()
  require_run_result("${_history_result}" "postprocess" "packaged ${_history_method}")
  if(_history_method STREQUAL "history_basis")
    string(JSON _history_basis GET "${_history_result}" archives 0)
    file(TO_CMAKE_PATH "${_history_basis}" _history_basis)
  endif()
endforeach()
