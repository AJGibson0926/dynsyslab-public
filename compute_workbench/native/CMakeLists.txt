cmake_minimum_required(VERSION 3.20)
project(dynsyslab_compute VERSION 0.1.0 LANGUAGES CXX)

if(MSVC)
  set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
endif()
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

include(CTest)
option(DYNSYSLAB_COMPUTE_BUILD_TESTS "Build native compute service tests" ON)
option(DYNSYSLAB_COMPUTE_BUILD_CLI "Build the standalone native workflow CLI" ON)

if(NOT TARGET dynsyslab_core)
  set(DYNSYSLAB_BUILD_PYTHON OFF CACHE BOOL "" FORCE)
  set(DYNSYSLAB_BUILD_CLI OFF CACHE BOOL "" FORCE)
  set(DYNSYSLAB_BUILD_TESTS OFF CACHE BOOL "" FORCE)
  set(DYNSYSLAB_ENABLE_ARCHIVE_IO ON CACHE BOOL "" FORCE)
  set(DYNSYSLAB_ENABLE_LINALG ON CACHE BOOL "" FORCE)
  set(DYNSYSLAB_ENABLE_SUNDIALS ON CACHE BOOL "" FORCE)
  set(DYNSYSLAB_ENABLE_SQLITE ON CACHE BOOL "" FORCE)
  set(DYNSYSLAB_BUILD_SHARED_MODULES ON CACHE BOOL "" FORCE)
  add_subdirectory(../../native ${CMAKE_CURRENT_BINARY_DIR}/core)
endif()
if(NOT TARGET dynsyslab_workflow)
  add_subdirectory(../../workflow/native ${CMAKE_CURRENT_BINARY_DIR}/workflow)
endif()
if(NOT TARGET dynsyslab_ui_support)
  add_subdirectory(../../native_ui ${CMAKE_CURRENT_BINARY_DIR}/ui-support)
endif()

if(DYNSYSLAB_COMPUTE_BUILD_CLI)
add_executable(dynsyslab_cli src/Cli.cpp)
  target_link_libraries(dynsyslab_cli PRIVATE dynsyslab_queue)
  target_compile_features(dynsyslab_cli PRIVATE cxx_std_17)
  if(MSVC)
    target_compile_options(dynsyslab_cli PRIVATE /W4 /permissive- /EHsc /utf-8)
  endif()
endif()
set_target_properties(dynsyslab_cli PROPERTIES OUTPUT_NAME "dynsyslab-cli")

add_executable(dynsyslab_compute WIN32
  src/main.cpp src/Registry.cpp src/Backend.cpp src/PanelCommon.cpp src/InlineValueControl.cpp
  src/SolverPanel.cpp src/SolverPanelLayout.cpp src/SolverPageChrome.cpp
  src/SolverPanelState.cpp src/SolverPanelSummary.cpp
  src/SolverPanelController.cpp src/SolverExecutionPane.cpp
  src/SolverExecutionContract.cpp
  src/ScientificParameterGrid.cpp src/ScientificParameterFormatting.cpp
  src/WaveformPreview.cpp src/SolveJobModel.cpp
  src/SolveJobPane.cpp src/SolveContextPanel.cpp src/WorkbenchSolveJobs.cpp
  src/SystemParameterDialog.cpp src/SettingsMerge.cpp
  src/KoopmanPipelineDialog.cpp src/CampaignDialog.cpp
  src/InitialConditionDialog.cpp src/ControllerSettingsDialog.cpp
  src/LegacyImportDialog.cpp
  src/KromControlSearchDialog.cpp
  src/LatticeInitialConditionDialog.cpp
  src/KromTrainingDialog.cpp
  src/QueueDialog.cpp src/ArchiveToolsDialog.cpp
  src/KromExtractionDialog.cpp
  src/ControlRegressionDialog.cpp
  src/KromValidationDialog.cpp
  src/BoundaryEigenfunctionDialog.cpp
  src/EigenvalueOptimizationDialog.cpp
  src/EstimationDialog.cpp
  src/DecompositionDialog.cpp
  src/DmdControlDialog.cpp
  src/PostprocessPanel.cpp src/Workbench.cpp src/WorkbenchArguments.cpp
  src/WorkbenchCheck.cpp src/WorkbenchChrome.cpp
  src/WorkbenchLayout.cpp)
if(WIN32)
  target_sources(dynsyslab_compute PRIVATE
    ../../native_ui/windows/DynSysLabAppIcon.rc)
endif()
set_target_properties(dynsyslab_compute PROPERTIES OUTPUT_NAME "DynSysLab")

target_include_directories(dynsyslab_compute PRIVATE include)
target_link_libraries(dynsyslab_compute PRIVATE dynsyslab_workflow dynsyslab_queue
  dynsyslab_ui_support)
if(WIN32)
  target_link_libraries(dynsyslab_compute PRIVATE comdlg32)
endif()
if(MSVC)
  target_compile_options(dynsyslab_compute PRIVATE /W4 /permissive- /EHsc /utf-8)
else()
  target_compile_options(dynsyslab_compute PRIVATE -Wall -Wextra -Wpedantic)
endif()

install(TARGETS dynsyslab_cli dynsyslab_compute RUNTIME DESTINATION bin)

if(BUILD_TESTING AND DYNSYSLAB_COMPUTE_BUILD_TESTS)
  add_executable(dynsyslab_native_solve_tests tests/native_solve_tests.cpp)
  target_link_libraries(dynsyslab_native_solve_tests PRIVATE dynsyslab_workflow)
  add_test(NAME dynsyslab_native_solve_tests COMMAND dynsyslab_native_solve_tests)
  add_executable(dynsyslab_krom_control_search_tests
    tests/krom_control_search_tests.cpp)
  target_link_libraries(dynsyslab_krom_control_search_tests PRIVATE dynsyslab_workflow)
  add_test(NAME dynsyslab_krom_control_search_tests
    COMMAND dynsyslab_krom_control_search_tests)
  set_tests_properties(dynsyslab_krom_control_search_tests PROPERTIES
    LABELS "native;workflow;ui;klqr;campaign" TIMEOUT 30)
  add_executable(dynsyslab_solve_job_model_tests
    tests/solve_job_model_tests.cpp src/SolveJobModel.cpp)
  target_include_directories(dynsyslab_solve_job_model_tests PRIVATE include)
  target_compile_features(dynsyslab_solve_job_model_tests PRIVATE cxx_std_17)
  add_test(NAME dynsyslab_solve_job_model_tests
    COMMAND dynsyslab_solve_job_model_tests)
  set_tests_properties(dynsyslab_solve_job_model_tests PROPERTIES
    LABELS "native;ui;solver;hierarchy")
  add_executable(dynsyslab_solve_job_pane_tests
    tests/solve_job_pane_tests.cpp src/SolveJobPane.cpp src/SolveJobModel.cpp
    src/PanelCommon.cpp)
  target_include_directories(dynsyslab_solve_job_pane_tests PRIVATE include)
  target_link_libraries(dynsyslab_solve_job_pane_tests PRIVATE dynsyslab_ui_support)
  target_compile_features(dynsyslab_solve_job_pane_tests PRIVATE cxx_std_17)
  add_test(NAME dynsyslab_solve_job_pane_tests
    COMMAND dynsyslab_solve_job_pane_tests)
  set_tests_properties(dynsyslab_solve_job_pane_tests PROPERTIES
    LABELS "native;ui;solver;hierarchy")
  add_executable(dynsyslab_solve_context_tests
    tests/solve_context_tests.cpp src/SolveContextPanel.cpp
    src/SolveJobModel.cpp src/PanelCommon.cpp)
  target_include_directories(dynsyslab_solve_context_tests PRIVATE include)
  target_link_libraries(dynsyslab_solve_context_tests PRIVATE
    dynsyslab_workflow dynsyslab_ui_support)
  target_compile_features(dynsyslab_solve_context_tests PRIVATE cxx_std_17)
  add_test(NAME dynsyslab_solve_context_tests
    COMMAND dynsyslab_solve_context_tests)
  set_tests_properties(dynsyslab_solve_context_tests PROPERTIES
    LABELS "native;ui;solver;initial-conditions")
  add_executable(dynsyslab_scientific_parameter_formatting_tests
    tests/scientific_parameter_formatting_tests.cpp
    src/ScientificParameterFormatting.cpp)
  target_include_directories(dynsyslab_scientific_parameter_formatting_tests
    PRIVATE include)
  target_compile_features(dynsyslab_scientific_parameter_formatting_tests
    PRIVATE cxx_std_17)
  add_test(NAME dynsyslab_scientific_parameter_formatting_tests
    COMMAND dynsyslab_scientific_parameter_formatting_tests)
  set_tests_properties(dynsyslab_scientific_parameter_formatting_tests PROPERTIES
    LABELS "native;ui;solver;parameters")
  add_executable(dynsyslab_solver_execution_contract_tests
    tests/solver_execution_contract_tests.cpp src/SolverExecutionContract.cpp)
  target_include_directories(dynsyslab_solver_execution_contract_tests PRIVATE include)
  target_link_libraries(dynsyslab_solver_execution_contract_tests PRIVATE
    dynsyslab_foundation)
  target_compile_features(dynsyslab_solver_execution_contract_tests PRIVATE cxx_std_17)
  add_test(NAME dynsyslab_solver_execution_contract_tests
    COMMAND dynsyslab_solver_execution_contract_tests)
  set_tests_properties(dynsyslab_solver_execution_contract_tests PROPERTIES
    LABELS "native;ui;solver;workflow")
  add_executable(dynsyslab_embedded_c_api_tests ../../native/tests/c_api_tests.cpp)
  target_link_libraries(dynsyslab_embedded_c_api_tests PRIVATE dynsyslab_c_api)
  add_test(NAME dynsyslab_embedded_c_api_tests COMMAND dynsyslab_embedded_c_api_tests)
  add_executable(dynsyslab_embedded_job_store_tests ../../native/tests/job_store_tests.cpp)
  target_link_libraries(dynsyslab_embedded_job_store_tests PRIVATE dynsyslab_operations)
  add_test(NAME dynsyslab_embedded_job_store_tests COMMAND dynsyslab_embedded_job_store_tests)
  add_executable(dynsyslab_workflow_c_api_tests
    ../../workflow/native/tests/workflow_c_api_tests.cpp)
  target_link_libraries(dynsyslab_workflow_c_api_tests PRIVATE dynsyslab_workflow_c_api)
  add_test(NAME dynsyslab_workflow_c_api_tests COMMAND dynsyslab_workflow_c_api_tests)
  add_executable(dynsyslab_control_factory_tests
    ../../workflow/native/tests/control_factory_tests.cpp)
  target_link_libraries(dynsyslab_control_factory_tests PRIVATE dynsyslab_workflow)
  add_test(NAME dynsyslab_control_factory_tests COMMAND dynsyslab_control_factory_tests)
  add_test(NAME dynsyslab_cli_queue_contract
    COMMAND ${CMAKE_COMMAND}
      -DCLI=$<TARGET_FILE:dynsyslab_cli>
      -DJOB=${CMAKE_SOURCE_DIR}/examples/native/duffing-free-solve.json
      -DROOT=${CMAKE_BINARY_DIR}/cli-queue-contract
      -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/CliQueueContract.cmake)
  add_test(NAME dynsyslab_cli_legacy_import_contract
    COMMAND ${CMAKE_COMMAND}
      -DCLI=$<TARGET_FILE:dynsyslab_cli>
      -DSOURCE=${CMAKE_SOURCE_DIR}/examples/native/legacy-dynbatch-example.csv
      -DROOT=${CMAKE_BINARY_DIR}/cli-legacy-import-contract
      -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/CliLegacyImportContract.cmake)
  add_test(NAME dynsyslab_cli_dmd_control_contract
    COMMAND ${CMAKE_COMMAND}
      -DCLI=$<TARGET_FILE:dynsyslab_cli>
      -DSOURCE=${CMAKE_SOURCE_DIR}/examples/native/legacy-dynbatch-example.csv
      -DROOT=${CMAKE_BINARY_DIR}/cli-dmdc-contract
      -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/CliDmdControlContract.cmake)
  add_test(NAME dynsyslab_cli_initial_condition_contract
    COMMAND ${CMAKE_COMMAND}
      -DCLI=$<TARGET_FILE:dynsyslab_cli>
      -DJOB=${CMAKE_SOURCE_DIR}/examples/native/duffing-free-solve.json
      -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/CliInitialConditionContract.cmake)
  add_library(dynsyslab_embedded_test_extension SHARED
    ../../native/tests/test_extension.cpp)
  target_include_directories(dynsyslab_embedded_test_extension PRIVATE ../../native/include)
  add_executable(dynsyslab_embedded_extension_tests
    ../../native/tests/extension_tests.cpp)
  target_link_libraries(dynsyslab_embedded_extension_tests PRIVATE dynsyslab_foundation)
  add_test(NAME dynsyslab_embedded_extension_tests COMMAND
    dynsyslab_embedded_extension_tests $<TARGET_FILE:dynsyslab_embedded_test_extension>)
endif()
