if(NOT STAGE_DIR)
  message(FATAL_ERROR "STAGE_DIR is required")
endif()

function(require_run_result output expected_kind label)
  string(JSON _schema GET "${output}" schema)
  string(JSON _kind GET "${output}" kind)
  string(JSON _archive_count LENGTH "${output}" archives)
  if(NOT _schema STREQUAL "dynsyslab.workflow-result/v1" OR
     NOT _kind STREQUAL "${expected_kind}" OR _archive_count LESS 1)
    message(FATAL_ERROR "${label} returned an invalid workflow result")
  endif()
endfunction()

set(_bin "${STAGE_DIR}/bin")
foreach(_required DynSysLab.exe DynSysLabViewer.exe dynsyslab-cli.exe
                  dynsyslab_foundation.dll dynsyslab_delay.dll dynsyslab_lattice.dll dynsyslab_pde.dll dynsyslab_linalg.dll
                  dynsyslab_eigenfunctions.dll dynsyslab_identification.dll dynsyslab_observation.dll
                  dynsyslab_workflow.dll dynsyslab_workflow_c_api.dll dynsyslab_phase_space.dll DynSysLabViewerCli.exe)
  if(NOT EXISTS "${_bin}/${_required}")
    message(FATAL_ERROR "packaged runtime is missing ${_required}")
  endif()
endforeach()
file(MAKE_DIRECTORY "${STAGE_DIR}/ui-smoke-archives" "${STAGE_DIR}/ui-smoke-cache")
foreach(_ui DynSysLab.exe DynSysLabViewer.exe)
  execute_process(COMMAND "${_bin}/${_ui}" --smoke-test
    --project-root "${STAGE_DIR}" --archive-root "${STAGE_DIR}/ui-smoke-archives"
    --cache "${STAGE_DIR}/ui-smoke-cache"
    RESULT_VARIABLE _ui_status ERROR_VARIABLE _ui_error TIMEOUT 20)
  if(NOT _ui_status EQUAL 0)
    message(FATAL_ERROR "packaged ${_ui} initialization gate failed: ${_ui_error}")
  endif()
endforeach()
file(GLOB_RECURSE _package_entries LIST_DIRECTORIES true "${STAGE_DIR}/*")
foreach(_entry IN LISTS _package_entries)
  get_filename_component(_entry_name "${_entry}" NAME)
  string(TOLOWER "${_entry_name}" _entry_name_lower)
  if(_entry_name_lower STREQUAL "__pycache__" OR
     _entry_name_lower MATCHES ".*\\.pyc$" OR
     _entry_name_lower MATCHES ".*\\.py$" OR
     _entry_name_lower MATCHES "^python.*\\.(dll|exe)$")
    message(FATAL_ERROR "packaged runtime unexpectedly contains Python artifact: ${_entry}")
  endif()
endforeach()
execute_process(COMMAND "${_bin}/dynsyslab-cli.exe" capabilities
  RESULT_VARIABLE _cap_status OUTPUT_VARIABLE _capabilities ERROR_VARIABLE _cap_error)
if(NOT _cap_status EQUAL 0 OR NOT _capabilities MATCHES "\"python_required\":false")
  message(FATAL_ERROR "packaged CLI capability gate failed: ${_cap_error}")
endif()
if(NOT _capabilities MATCHES "\"klqr\"" OR
   NOT _capabilities MATCHES "\"sdre\"")
  message(FATAL_ERROR "packaged CLI is missing native Koopman/classical control capabilities")
endif()
execute_process(COMMAND "${_bin}/dynsyslab-cli.exe" doctor
  --archive-root "${STAGE_DIR}/smoke-archives"
  RESULT_VARIABLE _doctor_status OUTPUT_VARIABLE _doctor ERROR_VARIABLE _doctor_error)
if(NOT _doctor_status EQUAL 0 OR NOT _doctor MATCHES "\"ok\":true")
  message(FATAL_ERROR "packaged CLI diagnostics gate failed: ${_doctor_error}")
endif()
set(_job "${STAGE_DIR}/smoke-job.json")
file(WRITE "${_job}" [=[{
  "schema":"dynsyslab.compute-job/v1","kind":"solve","label":"package-smoke",
  "system":{"name":"duffing","preset":"double_well_reference","parameters":{}},
  "time":{"start":0,"stop":0.01,"output_step":0.005,"coordinate_system":"automatic"},
  "initial_conditions":{"mode":"default","count":1,"clusters":1,"seed":0,"values":[],"ranges":[]},
  "controller":{"kind":"none","settings":{},"target":{"kind":"constant","state":[0,0]}},
  "solver":{"method":"RK45","rtol":1e-8,"atol":1e-10},
  "archive":{"enabled":true}
}]=])
execute_process(COMMAND "${_bin}/dynsyslab-cli.exe" run "${_job}"
  --archive-root "${STAGE_DIR}/smoke-archives"
  RESULT_VARIABLE _solve_status OUTPUT_VARIABLE _solve_output ERROR_VARIABLE _solve_error)
if(NOT _solve_status EQUAL 0)
  message(FATAL_ERROR "packaged solve gate failed: ${_solve_error}")
endif()
require_run_result("${_solve_output}" "solve" "packaged solve")
include("${CMAKE_CURRENT_LIST_DIR}/HistoryPackageSmoke.cmake")
string(JSON _solve_archive GET "${_solve_output}" archives 0)
file(TO_CMAKE_PATH "${_solve_archive}" _solve_archive_json)
set(_post_job "${STAGE_DIR}/postprocess-smoke-job.json")
file(WRITE "${_post_job}" "{\n"
  "  \"schema\":\"dynsyslab.compute-job/v1\",\"kind\":\"postprocess\",\n"
  "  \"label\":\"postprocess-package-smoke\",\"sources\":[\"${_solve_archive_json}\"],\n"
  "  \"methods\":[\n"
  "    {\"name\":\"bagged_pca\",\"enabled\":true,\"settings\":{\"rank\":1,\"bags\":3,\"seed\":7}},\n"
  "    {\"name\":\"robust_pca\",\"enabled\":true,\"settings\":{\"max_iterations\":5}}\n"
  "  ],\"archive\":{\"enabled\":true}\n"
  "}\n")
execute_process(COMMAND "${_bin}/dynsyslab-cli.exe" run "${_post_job}"
  --archive-root "${STAGE_DIR}/postprocess-smoke-archives"
  RESULT_VARIABLE _post_status OUTPUT_VARIABLE _post_output ERROR_VARIABLE _post_error)
if(NOT _post_status EQUAL 0)
  message(FATAL_ERROR "packaged post-processing gate failed: ${_post_error}")
endif()
require_run_result("${_post_output}" "postprocess" "packaged post-processing")
string(JSON _post_archive_count LENGTH "${_post_output}" archives)
if(NOT _post_archive_count EQUAL 2)
  message(FATAL_ERROR "packaged post-processing did not publish both derived archives")
endif()

set(_lqr_job "${STAGE_DIR}/lqr-smoke-job.json")
file(WRITE "${_lqr_job}" [=[{
  "schema":"dynsyslab.compute-job/v1","kind":"solve","label":"lqr-package-smoke",
  "system":{"name":"duffing","preset":"double_well_reference","parameters":{}},
  "time":{"start":0,"stop":0.01,"output_step":0.005,"coordinate_system":"automatic"},
  "initial_conditions":{"mode":"explicit","count":1,"clusters":1,"seed":0,
                        "values":[[-1,0.05]],"ranges":[]},
  "controller":{"kind":"lqr",
                "settings":{"Q":[10,1],"R":1,"linearization_state":[-1,0]},
                "target":{"kind":"constant","state":[-1,0]}},
  "solver":{"method":"RK4","step_size":0.001,
            "control_evaluation":"automatic"},
  "archive":{"enabled":true}
}]=])
execute_process(COMMAND "${_bin}/dynsyslab-cli.exe" run "${_lqr_job}"
  --archive-root "${STAGE_DIR}/lqr-smoke-archives"
  RESULT_VARIABLE _lqr_status OUTPUT_VARIABLE _lqr_output ERROR_VARIABLE _lqr_error)
if(NOT _lqr_status EQUAL 0)
  message(FATAL_ERROR "packaged LQR solve gate failed: ${_lqr_error}")
endif()
require_run_result("${_lqr_output}" "solve" "packaged LQR solve")

set(_delay_job "${STAGE_DIR}/delay-smoke-job.json")
file(WRITE "${_delay_job}" [=[{
  "schema":"dynsyslab.compute-job/v1","kind":"solve","label":"delay-package-smoke",
  "system":{"name":"mackey_glass","family":"dde","preset":"classic_chaos",
            "parameters":{"history_value":1.2}},
  "time":{"start":0,"stop":0.1,"output_step":0.05,"coordinate_system":"automatic"},
  "initial_conditions":{"mode":"explicit","count":2,"clusters":1,"seed":0,
                        "values":[[1.2],[1.0]],"ranges":[]},
  "controller":{"kind":"none","settings":{},"target":{"kind":"constant","state":[0]}},
  "solver":{"method":"method_of_steps","step_size":0.05},
  "archive":{"enabled":true,"data_role":"free_training"}
}]=])
execute_process(COMMAND "${_bin}/dynsyslab-cli.exe" run "${_delay_job}"
  --archive-root "${STAGE_DIR}/delay-smoke-archives"
  RESULT_VARIABLE _delay_status OUTPUT_VARIABLE _delay_output ERROR_VARIABLE _delay_error)
if(NOT _delay_status EQUAL 0)
  message(FATAL_ERROR "packaged delay solve gate failed: ${_delay_error}")
endif()
require_run_result("${_delay_output}" "solve" "packaged delay solve")

set(_lattice_job "${STAGE_DIR}/lattice-smoke-job.json")
file(WRITE "${_lattice_job}" [=[{
  "schema":"dynsyslab.compute-job/v1","kind":"solve","label":"lattice-package-smoke",
  "system":{"name":"potts","family":"lattice","preset":"q3_reference",
            "parameters":{"size":4,"states":3}},
  "time":{"start":0,"stop":3,"output_step":1,"coordinate_system":"automatic"},
  "initial_conditions":{"mode":"default","count":1,"clusters":1,"seed":11,
                        "values":[],"ranges":[]},
  "controller":{"kind":"none","settings":{},"target":{"kind":"constant","state":[0,0]}},
  "campaign":{"enabled":true,"kind":"temperature_sweep","mode":"cartesian",
              "axes":[{"path":"system.parameters.temperature","values":[1.0,2.0]}]},
  "solver":{"method":"metropolis"},
  "archive":{"enabled":true,"data_role":"free_training"}
}]=])
execute_process(COMMAND "${_bin}/dynsyslab-cli.exe" run "${_lattice_job}"
  --archive-root "${STAGE_DIR}/lattice-smoke-archives"
  RESULT_VARIABLE _lattice_status OUTPUT_VARIABLE _lattice_output
  ERROR_VARIABLE _lattice_error)
if(NOT _lattice_status EQUAL 0)
  message(FATAL_ERROR "packaged lattice campaign gate failed: ${_lattice_error}")
endif()
require_run_result("${_lattice_output}" "campaign" "packaged lattice campaign")

set(_pde_job "${STAGE_DIR}/pde-smoke-job.json")
file(WRITE "${_pde_job}" [=[{
  "schema":"dynsyslab.compute-job/v1","kind":"solve","label":"pde-package-smoke",
  "system":{"name":"heat_equation_1d","family":"pde",
            "preset":"periodic_decay","parameters":{"points":16,"diffusivity":0.001}},
  "time":{"start":0,"stop":0.02,"output_step":0.01,"coordinate_system":"automatic"},
  "initial_conditions":{"mode":"default","count":1,"clusters":1,"seed":0,
                        "values":[],"ranges":[]},
  "controller":{"kind":"none","settings":{},"target":{"kind":"constant","state":[]}},
  "solver":{"method":"RK4","step_size":0.005,
            "parallel":{"enabled":true,"partition_shape":[4],
                         "ghost_width":1,"local_workers":2}},
  "archive":{"enabled":true,"data_role":"free_training"}
}]=])
execute_process(COMMAND "${_bin}/dynsyslab-cli.exe" run "${_pde_job}"
  --archive-root "${STAGE_DIR}/pde-smoke-archives"
  RESULT_VARIABLE _pde_status OUTPUT_VARIABLE _pde_output ERROR_VARIABLE _pde_error)
if(NOT _pde_status EQUAL 0)
  message(FATAL_ERROR "packaged PDE solve gate failed: ${_pde_error}")
endif()
require_run_result("${_pde_output}" "solve" "packaged PDE solve")

set(_dedup_root "${STAGE_DIR}/deduplication-smoke-archives")
file(MAKE_DIRECTORY "${_dedup_root}")
execute_process(COMMAND "${_bin}/dynsyslab-cli.exe" archive deduplicate
  --archive-root "${_dedup_root}" --all-ages
  RESULT_VARIABLE _dedup_status OUTPUT_VARIABLE _dedup_output ERROR_VARIABLE _dedup_error)
if(NOT _dedup_status EQUAL 0)
  message(FATAL_ERROR "packaged deduplication gate failed: ${_dedup_error}")
endif()
string(JSON _dedup_schema GET "${_dedup_output}" schema)
string(JSON _dedup_dry_run GET "${_dedup_output}" dry_run)
string(JSON _dedup_recoverable GET "${_dedup_output}" recoverable)
if(NOT _dedup_schema STREQUAL "dynsyslab.archive-deduplication/v1" OR
   NOT _dedup_dry_run OR NOT _dedup_recoverable)
  message(FATAL_ERROR "packaged deduplication safety contract failed")
endif()

set(_sdk_source "${STAGE_DIR}/resources/sdk/examples/native-extension")
set(_sdk_build "${STAGE_DIR}/sdk-extension-build")
file(REMOVE_RECURSE "${_sdk_build}")
execute_process(COMMAND "${CMAKE_COMMAND}" -S "${_sdk_source}" -B "${_sdk_build}"
  -G Ninja -DDYNSYSLAB_SDK_DIR=${STAGE_DIR}
  RESULT_VARIABLE _sdk_configure_status OUTPUT_VARIABLE _sdk_configure_output
  ERROR_VARIABLE _sdk_configure_error)
if(NOT _sdk_configure_status EQUAL 0)
  message(FATAL_ERROR "installed extension SDK configure failed: ${_sdk_configure_error}")
endif()
execute_process(COMMAND "${CMAKE_COMMAND}" --build "${_sdk_build}"
  RESULT_VARIABLE _sdk_build_status OUTPUT_VARIABLE _sdk_build_output
  ERROR_VARIABLE _sdk_build_error)
if(NOT _sdk_build_status EQUAL 0)
  message(FATAL_ERROR "installed extension SDK build failed: ${_sdk_build_error}")
endif()
set(_extension "${_sdk_build}/dynsyslab-example-oscillator.dll")
if(NOT EXISTS "${_extension}")
  message(FATAL_ERROR "installed extension SDK did not produce the example DLL")
endif()
file(TO_CMAKE_PATH "${_extension}" _extension_json)
set(_extension_job "${STAGE_DIR}/extension-smoke-job.json")
file(WRITE "${_extension_job}" "{\n"
  "  \"schema\":\"dynsyslab.compute-job/v1\",\"kind\":\"solve\",\n"
  "  \"label\":\"extension-smoke\",\"extensions\":[\"${_extension_json}\"],\n"
  "  \"system\":{\"name\":\"example_oscillator\",\"parameters\":{}},\n"
  "  \"time\":{\"start\":0,\"stop\":0.01,\"output_step\":0.005},\n"
  "  \"initial_conditions\":{\"mode\":\"default\",\"count\":1,\"clusters\":1,\"values\":[],\"ranges\":[]},\n"
  "  \"controller\":{\"kind\":\"none\",\"settings\":{},\"target\":{\"kind\":\"constant\",\"state\":[0,0]}},\n"
  "  \"solver\":{\"method\":\"RK4\"},\"archive\":{\"enabled\":true}\n"
  "}\n")
execute_process(COMMAND "${_bin}/dynsyslab-cli.exe" run "${_extension_job}"
  --archive-root "${STAGE_DIR}/extension-smoke-archives"
  RESULT_VARIABLE _extension_status OUTPUT_VARIABLE _extension_output
  ERROR_VARIABLE _extension_error)
if(NOT _extension_status EQUAL 0)
  message(FATAL_ERROR "packaged extension solve gate failed: ${_extension_error}")
endif()
require_run_result("${_extension_output}" "solve" "packaged extension solve")
message(STATUS "DynSysLab isolated package smoke gate passed")

set(_phase_query "${STAGE_DIR}/phase-smoke-query.json")
file(WRITE "${_phase_query}" [=[{"schema":"dynsyslab.phase-space-query/v1","mode":"compute","system":{"name":"duffing"},"options":{"lower":[-2,-2],"upper":[2,2],"periodic_orbits":false,"manifolds":false}}]=])
execute_process(COMMAND "${_bin}/dynsyslab-cli.exe" phase-space "${_phase_query}" --archive-root "${STAGE_DIR}/phase-smoke-archives"
    RESULT_VARIABLE _phase_status OUTPUT_VARIABLE _phase_output ERROR_VARIABLE _phase_error TIMEOUT 60)
if(NOT _phase_status EQUAL 0 OR NOT _phase_output MATCHES "dynsyslab.phase-space/v1")
    message(FATAL_ERROR "Packaged phase-space API failed: ${_phase_error}")
endif()
set(_research_query "${STAGE_DIR}/research-smoke-query.json")
file(WRITE "${_research_query}" [=[{"schema":"dynsyslab.research-query/v1","action":"list","system":{"name":"rpe"}}]=])
execute_process(COMMAND "${_bin}/dynsyslab-cli.exe" research "${_research_query}" --archive-root "${STAGE_DIR}/phase-smoke-archives"
    RESULT_VARIABLE _research_status OUTPUT_VARIABLE _research_output ERROR_VARIABLE _research_error TIMEOUT 30)
if(NOT _research_status EQUAL 0 OR NOT _research_output MATCHES "chang-chen-1986")
    message(FATAL_ERROR "Packaged system bibliography API failed: ${_research_error}")
endif()
message(STATUS "Packaged phase-space and system bibliography APIs passed")

# Analytic factories belong to non-ODE systems and must be available from the packaged CLI.
set(_heat_query "${STAGE_DIR}/analytic-heat-query.json")
file(WRITE "${_heat_query}" [=[{"schema":"dynsyslab.analytic-query/v1","system":{"name":"heat_equation_1d","parameters":{"points":3}},"initial_state":[0,1,-1],"time":[0,0.1]}]=])
set(_delay_query "${STAGE_DIR}/analytic-delay-query.json")
file(WRITE "${_delay_query}" [=[{"schema":"dynsyslab.analytic-query/v1","system":{"name":"mackey_glass","parameters":{"exponent":0,"beta":0.3,"delay":1}},"initial_state":[1],"time":[0,2.5]}]=])
foreach(_analytic_query IN ITEMS "${_heat_query}" "${_delay_query}")
  execute_process(COMMAND "${_bin}/dynsyslab-cli.exe" analytic "${_analytic_query}"
    RESULT_VARIABLE _analytic_status OUTPUT_VARIABLE _analytic_output ERROR_VARIABLE _analytic_error TIMEOUT 20)
  if(NOT _analytic_status EQUAL 0)
    message(FATAL_ERROR "packaged non-ODE analytic query failed: ${_analytic_error}")
  endif()
  string(JSON _analytic_complete GET "${_analytic_output}" status)
  string(JSON _analytic_count LENGTH "${_analytic_output}" values)
  if(NOT _analytic_complete STREQUAL "complete" OR _analytic_count LESS 2)
    message(FATAL_ERROR "packaged non-ODE analytic query produced no complete reference")
  endif()
endforeach()
message(STATUS "Packaged heat and delay analytic families passed")

set(_mechanical_query "${STAGE_DIR}/analytic-mechanical-query.json")
file(WRITE "${_mechanical_query}" [=[{"schema":"dynsyslab.analytic-query/v1","system":{"name":"double_pendulum","parameters":{"gravity_m_s2":0,"damping1_n_m_s":0,"damping2_n_m_s":0}},"initial_state":[0.6,0.5,0,0.5],"time":[0,1,8]}]=])
execute_process(COMMAND "${_bin}/dynsyslab-cli.exe" analytic "${_mechanical_query}"
  RESULT_VARIABLE _mechanical_status OUTPUT_VARIABLE _mechanical_output ERROR_VARIABLE _mechanical_error TIMEOUT 20)
if(NOT _mechanical_status EQUAL 0)
  message(FATAL_ERROR "Packaged double-pendulum analytic family failed: ${_mechanical_error}")
endif()
string(JSON _mechanical_complete GET "${_mechanical_output}" status)
if(NOT _mechanical_complete STREQUAL "complete")
  message(FATAL_ERROR "Packaged double-pendulum reference incomplete")
endif()
set(_lauricella_query "${STAGE_DIR}/lauricella-query.json")
file(WRITE "${_lauricella_query}" [=[{"schema":"dynsyslab.special-function-query/v1","function":"lauricella_fd","a":1,"b":[1],"c":2,"x":[0.5]}]=])
execute_process(COMMAND "${_bin}/dynsyslab-cli.exe" special-function "${_lauricella_query}"
  RESULT_VARIABLE _lauricella_status OUTPUT_VARIABLE _lauricella_output ERROR_VARIABLE _lauricella_error TIMEOUT 20)
if(NOT _lauricella_status EQUAL 0 OR NOT _lauricella_output MATCHES "1.386294")
  message(FATAL_ERROR "Packaged Lauricella evaluator failed: ${_lauricella_error}")
endif()
message(STATUS "Packaged double-pendulum and Lauricella APIs passed")


# Catalogue limits must survive the isolated package, including custom bubble classes.
foreach(_catalogue_case IN ITEMS kuramoto-pair hh-constant-voltage spring-pendulum-synchronous marmottant-stationary reddy-szeri-shapes)
  execute_process(COMMAND "${_bin}/dynsyslab-cli.exe" analytic "${CMAKE_CURRENT_LIST_DIR}/../examples/native/${_catalogue_case}-analytic-query.json"
    RESULT_VARIABLE _catalogue_status OUTPUT_VARIABLE _catalogue_output ERROR_VARIABLE _catalogue_error TIMEOUT 20)
  if(NOT _catalogue_status EQUAL 0)
    message(FATAL_ERROR "Packaged ${_catalogue_case} query failed: ${_catalogue_error}")
  endif()
  string(JSON _catalogue_complete GET "${_catalogue_output}" status)
  if(NOT _catalogue_complete STREQUAL "complete")
    message(FATAL_ERROR "Packaged ${_catalogue_case} reference incomplete")
  endif()
endforeach()
message(STATUS "Packaged catalogue analytic limits passed")


execute_process(COMMAND "${_bin}/dynsyslab-cli.exe" phase-space "${CMAKE_CURRENT_LIST_DIR}/../examples/native/rpe-periodic-phase-query.json"
  --archive-root "${STAGE_DIR}/periodic-phase-smoke-archives"
  RESULT_VARIABLE _forced_status OUTPUT_VARIABLE _forced_output ERROR_VARIABLE _forced_error TIMEOUT 30)
if(NOT _forced_status EQUAL 0)
  message(FATAL_ERROR "Packaged forced RPE sidecar failed: ${_forced_error}")
endif()
string(JSON _forced_count LENGTH "${_forced_output}" periodic_orbits)
string(JSON _forced_mode GET "${_forced_output}" options dynamics)
if(_forced_count LESS 1 OR NOT _forced_mode STREQUAL "periodic-forcing")
  message(FATAL_ERROR "Packaged forced RPE sidecar omitted the periodic response")
endif()
message(STATUS "Packaged periodic-forcing sidecar passed")
