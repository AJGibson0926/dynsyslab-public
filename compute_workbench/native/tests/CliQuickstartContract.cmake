# dynsyslab-cli quickstart: welcome, bundled examples, results table, figures and next steps
# as Markdown; distinct archives on repeated runs; JSON form; nothing written to the working
# directory.
file(REMOVE_RECURSE "${ROOT}")
file(MAKE_DIRECTORY "${ROOT}/data/archives" "${ROOT}/cwd")

execute_process(COMMAND "${CLI}" quickstart --archive-root "${ROOT}/data/archives"
  WORKING_DIRECTORY "${ROOT}/cwd" RESULT_VARIABLE code OUTPUT_VARIABLE markdown ERROR_VARIABLE error)
if(NOT code EQUAL 0)
  message(FATAL_ERROR "quickstart failed: ${error}")
endif()
foreach(expected "# Welcome to DynSysLab" "## Quickstart results" "| Duffing double well |"
                 "| Lorenz 1963 (chaotic) |" "## Figures" "## Things you can ask next"
                 "The desktop applications were not started")
  string(FIND "${markdown}" "${expected}" at)
  if(at EQUAL -1)
    message(FATAL_ERROR "quickstart Markdown lacks '${expected}':\n${markdown}")
  endif()
endforeach()
string(FIND "${markdown}" "Agents:" agents)
if(NOT agents EQUAL -1)
  message(FATAL_ERROR "quickstart output addresses agents directly")
endif()

execute_process(COMMAND "${CLI}" quickstart --archive-root "${ROOT}/data/archives" --format json
  WORKING_DIRECTORY "${ROOT}/cwd" RESULT_VARIABLE code OUTPUT_VARIABLE report ERROR_VARIABLE error)
if(NOT code EQUAL 0)
  message(FATAL_ERROR "quickstart --format json failed: ${error}")
endif()
string(JSON schema GET "${report}" schema)
string(JSON runs LENGTH "${report}" runs)
if(NOT schema STREQUAL "dynsyslab.quickstart/v1" OR runs LESS 2)
  message(FATAL_ERROR "quickstart JSON report is incomplete: ${report}")
endif()

file(GLOB duffing "${ROOT}/data/archives/solutions/duffing/quickstart-duffing-double-well__*.npz")
list(LENGTH duffing count)
if(NOT count EQUAL 2)
  message(FATAL_ERROR "two quickstart runs should leave two distinct Duffing archives, found ${count}")
endif()
file(GLOB written "${ROOT}/cwd/*")
if(written)
  message(FATAL_ERROR "quickstart wrote into the working directory: ${written}")
endif()
