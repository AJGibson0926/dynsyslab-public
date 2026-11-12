if(NOT EXISTS "${CLI}")
  message(FATAL_ERROR "native CLI is missing")
endif()
if(NOT EXISTS "${SOURCE}")
  message(FATAL_ERROR "legacy CSV fixture is missing")
endif()
file(REMOVE_RECURSE "${ROOT}")
execute_process(
  COMMAND "${CLI}" import legacy-csv "${SOURCE}"
    --archive-root "${ROOT}" --map-system duffing
    --preset double_well_reference --coordinates dimensionless
    --label cli-legacy-import
  RESULT_VARIABLE status
  OUTPUT_VARIABLE output
  ERROR_VARIABLE error)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "legacy import failed: ${error}")
endif()
string(FIND "${output}" "dynsyslab.legacy-csv-import-result/v1" schema)
if(schema EQUAL -1)
  message(FATAL_ERROR "legacy import result schema is missing: ${output}")
endif()
file(GLOB_RECURSE imported "${ROOT}/solutions/*.npz")
list(LENGTH imported imported_count)
if(NOT imported_count EQUAL 1)
  message(FATAL_ERROR "legacy import did not publish exactly one trajectory")
endif()
file(REMOVE_RECURSE "${ROOT}")
