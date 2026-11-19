cmake_minimum_required(VERSION 3.20)
project(dynsyslab_viewer VERSION 0.1.0 LANGUAGES CXX)

if(MSVC)
  set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
endif()
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

include(CTest)
option(DYNSYSLAB_VIEWER_BUILD_TESTS "Build native Viewer service tests" ON)

if(NOT TARGET dynsyslab_core)
  set(DYNSYSLAB_BUILD_PYTHON OFF CACHE BOOL "" FORCE)
  set(DYNSYSLAB_BUILD_CLI OFF CACHE BOOL "" FORCE)
  set(DYNSYSLAB_BUILD_TESTS OFF CACHE BOOL "" FORCE)
  set(DYNSYSLAB_ENABLE_ARCHIVE_IO ON CACHE BOOL "" FORCE)
  # Standalone Viewer only decodes archives; it does not compile compute kernels.
  set(DYNSYSLAB_ENABLE_LINALG OFF CACHE BOOL "" FORCE)
  set(DYNSYSLAB_ENABLE_SUNDIALS OFF CACHE BOOL "" FORCE)
  set(DYNSYSLAB_ENABLE_SQLITE ON CACHE BOOL "" FORCE)
  set(DYNSYSLAB_BUILD_SHARED_MODULES ON CACHE BOOL "" FORCE)
  add_subdirectory(../../native ${CMAKE_CURRENT_BINARY_DIR}/core)
endif()
if(NOT TARGET dynsyslab_ui_support)
  add_subdirectory(../../native_ui ${CMAKE_CURRENT_BINARY_DIR}/ui-support)
endif()

add_library(dynsyslab_viewer_services STATIC
  src/Backend.cpp
  src/ArchiveSystem.cpp
  src/PhaseContextBackend.cpp
  src/PhaseSidecarBridge.cpp
  src/PhaseContextOverlay.cpp
  src/PhaseSurfaceOverlay.cpp
  src/ViewerSession.cpp
  src/ViewerWindow.cpp
  src/PlotSnapshot.cpp
  src/PlotRender.cpp
  src/TimeSeriesOverlay.cpp
  src/RecordWindow.cpp
  src/ComparisonPayload.cpp
  src/DisplayCoordinates.cpp
  src/TrajectoryPresentation.cpp
  src/AnalyticOverlay.cpp
  src/FieldPayload.cpp
  src/FourierPayload.cpp
  src/HierarchyPayload.cpp
  src/LatticePayload.cpp
  src/NativePayload.cpp
  src/ArchiveMatrixPayload.cpp
  src/ObservablePayload.cpp
  src/ObservableWorkspace.cpp
  src/PanelCatalog.cpp
  src/KoopmanValidationPayload.cpp
  src/SindyDynamicsPayload.cpp
  src/EigenvalueOptimizationPayload.cpp
  src/SpectralMeasurePayload.cpp
  src/ComplexStemPlot.cpp
  src/ModeWeightPayload.cpp
  src/LinearSystemsPayload.cpp
  src/LinearSystemsLoopPayload.cpp
  src/KromControlCampaignPayload.cpp
  src/KromPayload.cpp
  src/KromTraining.cpp
  src/KromWorkspace.cpp)
target_include_directories(dynsyslab_viewer_services PUBLIC include)
target_link_libraries(dynsyslab_viewer_services PUBLIC
  dynsyslab_ui_support dynsyslab_core dynsyslab_koopman dynsyslab_operations)
target_compile_features(dynsyslab_viewer_services PUBLIC cxx_std_17)

add_executable(dynsyslab_viewer WIN32
  src/main.cpp
  src/AxisRangeDialog.cpp
  src/ComparisonDialog.cpp
  src/PanelPickerDialog.cpp
  src/PlotExport.cpp
  src/Workbench.cpp
  src/WorkbenchArguments.cpp
  src/WorkbenchCheck.cpp
  src/WorkbenchKromCheck.cpp
  src/WorkbenchCatalog.cpp
  src/WorkbenchWorkspace.cpp
  src/WorkbenchTasks.cpp
  src/WorkbenchPhaseContext.cpp
  src/WorkbenchSnapshot.cpp
  src/WorkbenchLayout.cpp)
target_sources(dynsyslab_viewer PRIVATE src/WorkbenchActions.cpp)
if(WIN32)
  target_sources(dynsyslab_viewer PRIVATE
    ../../native_ui/windows/DynSysLabViewerIcon.rc)
endif()
set_target_properties(dynsyslab_viewer PROPERTIES OUTPUT_NAME "DynSysLabViewer")

target_link_libraries(dynsyslab_viewer PRIVATE dynsyslab_viewer_services)
target_link_libraries(dynsyslab_viewer PRIVATE comdlg32 ole32)

if(MSVC)
  target_compile_options(dynsyslab_viewer PRIVATE /W4 /permissive- /EHsc /utf-8)
else()
  target_compile_options(dynsyslab_viewer PRIVATE -Wall -Wextra -Wpedantic)
endif()

install(TARGETS dynsyslab_viewer RUNTIME DESTINATION bin)

add_executable(dynsyslab_viewer_cli tools/ViewerCli.cpp tools/ViewerSessionCli.cpp)
set_target_properties(dynsyslab_viewer_cli PROPERTIES OUTPUT_NAME "DynSysLabViewerCli")
target_link_libraries(dynsyslab_viewer_cli PRIVATE dynsyslab_viewer_services)
target_compile_features(dynsyslab_viewer_cli PRIVATE cxx_std_17)
if(MSVC)
  target_compile_options(dynsyslab_viewer_cli PRIVATE /W4 /permissive- /EHsc /utf-8)
endif()
# The sidecar is loaded dynamically; keep full-suite builds/package installs complete.
if(TARGET dynsyslab_workflow_c_api)
  add_dependencies(dynsyslab_viewer dynsyslab_workflow_c_api)
  add_dependencies(dynsyslab_viewer_cli dynsyslab_workflow_c_api)
endif()
install(TARGETS dynsyslab_viewer_cli RUNTIME DESTINATION bin)
install(FILES include/dynsyslab_viewer/ViewerSession.hpp
  include/dynsyslab_viewer/TimeSeriesOverlay.hpp
  include/dynsyslab_viewer/KromWorkspace.hpp
  include/dynsyslab_viewer/ViewerWindow.hpp
  include/dynsyslab_viewer/Backend.hpp include/dynsyslab_viewer/PanelCatalog.hpp
  DESTINATION include/dynsyslab_viewer)

if(BUILD_TESTING AND DYNSYSLAB_VIEWER_BUILD_TESTS)
  add_executable(dynsyslab_time_series_overlay_tests tests/time_series_overlay_tests.cpp)
  target_link_libraries(dynsyslab_time_series_overlay_tests PRIVATE dynsyslab_viewer_services)
  add_test(NAME dynsyslab_time_series_overlay_tests COMMAND dynsyslab_time_series_overlay_tests)
  add_test(NAME dynsyslab_viewer_workbench_background_contract COMMAND ${CMAKE_COMMAND}
    -DVIEWER=$<TARGET_FILE:dynsyslab_viewer>
    -DVIEWER_CLI=$<TARGET_FILE:dynsyslab_viewer_cli>
    -DOBSERVABLE_FIXTURE=$<TARGET_FILE:dynsyslab_observable_workspace_tests>
    -DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}/workbench-background-test
    -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/viewer_workbench_background.cmake)
  set_tests_properties(dynsyslab_viewer_workbench_background_contract PROPERTIES
    LABELS "native;ui;headless;focus" TIMEOUT 90)
  add_test(NAME dynsyslab_viewer_cli_background_contract COMMAND ${CMAKE_COMMAND}
    -DVIEWER_CLI=$<TARGET_FILE:dynsyslab_viewer_cli>
    -DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}/cli-background-test
    -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/viewer_cli_background.cmake)
  set_tests_properties(dynsyslab_viewer_cli_background_contract PROPERTIES
    LABELS "native;ui;cli;headless")
  if(TARGET dynsyslab_workflow)
    add_executable(dynsyslab_phase_forcing_workflow_tests tests/phase_forcing_workflow_tests.cpp)
    target_link_libraries(dynsyslab_phase_forcing_workflow_tests PRIVATE dynsyslab_viewer_services dynsyslab_workflow)
    if(TARGET dynsyslab_workflow_c_api)
      add_dependencies(dynsyslab_phase_forcing_workflow_tests dynsyslab_workflow_c_api)
    endif()
    add_test(NAME dynsyslab_phase_forcing_workflow_tests COMMAND dynsyslab_phase_forcing_workflow_tests)
    add_executable(dynsyslab_analytic_catalogue_workflow_tests tests/analytic_catalogue_workflow_tests.cpp)
    target_link_libraries(dynsyslab_analytic_catalogue_workflow_tests PRIVATE dynsyslab_viewer_services dynsyslab_workflow)
    add_test(NAME dynsyslab_analytic_catalogue_workflow_tests COMMAND dynsyslab_analytic_catalogue_workflow_tests)
    add_executable(dynsyslab_analytic_mechanical_workflow_tests tests/analytic_mechanical_workflow_tests.cpp)
    target_link_libraries(dynsyslab_analytic_mechanical_workflow_tests PRIVATE dynsyslab_viewer_services dynsyslab_workflow)
    add_test(NAME dynsyslab_analytic_mechanical_workflow_tests COMMAND dynsyslab_analytic_mechanical_workflow_tests)
    add_executable(dynsyslab_nonode_analytic_workflow_tests tests/nonode_analytic_workflow_tests.cpp)
    target_link_libraries(dynsyslab_nonode_analytic_workflow_tests PRIVATE dynsyslab_viewer_services dynsyslab_workflow)
    add_test(NAME dynsyslab_nonode_analytic_workflow_tests COMMAND dynsyslab_nonode_analytic_workflow_tests)
  endif()
  add_executable(dynsyslab_phase_context_tests tests/phase_context_tests.cpp)
  target_link_libraries(dynsyslab_phase_context_tests PRIVATE dynsyslab_viewer_services)
  if(TARGET dynsyslab_workflow_c_api)
    add_dependencies(dynsyslab_phase_context_tests dynsyslab_workflow_c_api)
    add_test(NAME dynsyslab_phase_context_tests COMMAND dynsyslab_phase_context_tests)
  endif()
  add_executable(dynsyslab_plot_scale_scroll_tests tests/plot_scale_scroll_tests.cpp)
  target_link_libraries(dynsyslab_plot_scale_scroll_tests PRIVATE dynsyslab_viewer_services)
  add_test(NAME dynsyslab_plot_scale_scroll_tests COMMAND dynsyslab_plot_scale_scroll_tests)
  add_executable(dynsyslab_krom_workspace_tests tests/krom_workspace_tests.cpp)
  target_link_libraries(dynsyslab_krom_workspace_tests PRIVATE dynsyslab_viewer_services)
  add_test(NAME dynsyslab_krom_workspace_tests COMMAND dynsyslab_krom_workspace_tests)
  add_executable(dynsyslab_observable_workspace_tests tests/observable_workspace_tests.cpp)
  target_link_libraries(dynsyslab_observable_workspace_tests PRIVATE dynsyslab_viewer_services)
  add_test(NAME dynsyslab_observable_workspace_tests COMMAND dynsyslab_observable_workspace_tests)
  add_executable(dynsyslab_plot_presentation_tests tests/plot_presentation_tests.cpp)
  target_link_libraries(dynsyslab_plot_presentation_tests PRIVATE dynsyslab_viewer_services)
  add_test(NAME dynsyslab_plot_presentation_tests COMMAND dynsyslab_plot_presentation_tests)
  add_executable(dynsyslab_native_payload_tests tests/native_payload_tests.cpp)
  add_executable(dynsyslab_sindy_dynamics_payload_tests tests/sindy_dynamics_payload_tests.cpp)
  target_link_libraries(dynsyslab_sindy_dynamics_payload_tests PRIVATE dynsyslab_viewer_services)
  add_test(NAME dynsyslab_sindy_dynamics_payload_tests COMMAND dynsyslab_sindy_dynamics_payload_tests)
  target_link_libraries(dynsyslab_native_payload_tests PRIVATE dynsyslab_viewer_services)
  add_test(NAME dynsyslab_native_payload_tests COMMAND dynsyslab_native_payload_tests)
  add_executable(dynsyslab_linear_systems_payload_tests tests/linear_systems_payload_tests.cpp)
  target_link_libraries(dynsyslab_linear_systems_payload_tests PRIVATE dynsyslab_viewer_services)
  add_test(NAME dynsyslab_linear_systems_payload_tests COMMAND dynsyslab_linear_systems_payload_tests)
  add_executable(dynsyslab_viewer_display_policy_tests tests/display_policy_tests.cpp)
  target_link_libraries(dynsyslab_viewer_display_policy_tests PRIVATE dynsyslab_viewer_services)
  add_test(NAME dynsyslab_viewer_display_policy_tests COMMAND dynsyslab_viewer_display_policy_tests)
  add_executable(dynsyslab_hierarchy_catalog_tests tests/hierarchy_catalog_tests.cpp)
  target_link_libraries(dynsyslab_hierarchy_catalog_tests PRIVATE dynsyslab_viewer_services)
  add_test(NAME dynsyslab_hierarchy_catalog_tests COMMAND dynsyslab_hierarchy_catalog_tests)
  add_executable(dynsyslab_archive_focus_tests tests/archive_focus_tests.cpp)
  target_link_libraries(dynsyslab_archive_focus_tests PRIVATE dynsyslab_ui_support)
  add_test(NAME dynsyslab_archive_focus_tests COMMAND dynsyslab_archive_focus_tests)
  set_tests_properties(dynsyslab_archive_focus_tests PROPERTIES
    LABELS "native;ui;hierarchy;focus")
  add_executable(dynsyslab_catalog_presentation_tests
    tests/catalog_presentation_tests.cpp)
  target_link_libraries(dynsyslab_catalog_presentation_tests PRIVATE
    dynsyslab_ui_support)
  add_test(NAME dynsyslab_catalog_presentation_tests
    COMMAND dynsyslab_catalog_presentation_tests)
  set_tests_properties(dynsyslab_catalog_presentation_tests PROPERTIES
    LABELS "native;ui;hierarchy;presentation")
  add_executable(dynsyslab_catalog_tree_sort_tests
    tests/catalog_tree_sort_tests.cpp)
  target_link_libraries(dynsyslab_catalog_tree_sort_tests PRIVATE
    dynsyslab_ui_support)
  add_test(NAME dynsyslab_catalog_tree_sort_tests
    COMMAND dynsyslab_catalog_tree_sort_tests)
  set_tests_properties(dynsyslab_catalog_tree_sort_tests PROPERTIES
    LABELS "native;ui;hierarchy;sorting")
  add_executable(dynsyslab_runtime_paths_tests tests/runtime_paths_tests.cpp)
  target_link_libraries(dynsyslab_runtime_paths_tests PRIVATE dynsyslab_ui_support)
  add_test(NAME dynsyslab_runtime_paths_tests COMMAND dynsyslab_runtime_paths_tests)
  add_executable(dynsyslab_horizontal_pane_layout_tests
    tests/horizontal_pane_layout_tests.cpp)
  target_link_libraries(dynsyslab_horizontal_pane_layout_tests
    PRIVATE dynsyslab_ui_support)
  add_test(NAME dynsyslab_horizontal_pane_layout_tests
    COMMAND dynsyslab_horizontal_pane_layout_tests)
  add_executable(dynsyslab_ornamental_scrollbar_tests
    tests/ornamental_scrollbar_tests.cpp)
  target_link_libraries(dynsyslab_ornamental_scrollbar_tests PRIVATE
    dynsyslab_ui_support)
  add_test(NAME dynsyslab_ornamental_scrollbar_tests
    COMMAND dynsyslab_ornamental_scrollbar_tests)
  set_tests_properties(dynsyslab_ornamental_scrollbar_tests PROPERTIES
    LABELS "native;ui;scrolling")
  add_executable(dynsyslab_ornamental_combo_box_tests
    tests/ornamental_combo_box_tests.cpp)
  target_link_libraries(dynsyslab_ornamental_combo_box_tests PRIVATE
    dynsyslab_ui_support)
  add_test(NAME dynsyslab_ornamental_combo_box_tests
    COMMAND dynsyslab_ornamental_combo_box_tests)
  set_tests_properties(dynsyslab_ornamental_combo_box_tests PROPERTIES
    LABELS "native;ui;controls")
  add_executable(dynsyslab_koopman_validation_payload_tests
    tests/koopman_validation_payload_tests.cpp)
  target_link_libraries(dynsyslab_koopman_validation_payload_tests
    PRIVATE dynsyslab_viewer_services)
  add_test(NAME dynsyslab_koopman_validation_payload_tests
    COMMAND dynsyslab_koopman_validation_payload_tests)
  add_executable(dynsyslab_eigenvalue_optimization_payload_tests
    tests/eigenvalue_optimization_payload_tests.cpp)
  target_link_libraries(dynsyslab_eigenvalue_optimization_payload_tests
    PRIVATE dynsyslab_viewer_services)
  add_test(NAME dynsyslab_eigenvalue_optimization_payload_tests
    COMMAND dynsyslab_eigenvalue_optimization_payload_tests)
  add_executable(dynsyslab_krom_control_campaign_payload_tests
    tests/krom_control_campaign_payload_tests.cpp)
  target_link_libraries(dynsyslab_krom_control_campaign_payload_tests
    PRIVATE dynsyslab_viewer_services)
  add_test(NAME dynsyslab_krom_control_campaign_payload_tests
    COMMAND dynsyslab_krom_control_campaign_payload_tests)
  add_executable(dynsyslab_plot_export_tests
    tests/plot_export_tests.cpp src/PlotExport.cpp)
  target_link_libraries(dynsyslab_plot_export_tests PRIVATE
    dynsyslab_viewer_services comdlg32 ole32)
  add_test(NAME dynsyslab_plot_export_tests COMMAND dynsyslab_plot_export_tests)
  add_executable(dynsyslab_panel_workspace_tests tests/panel_workspace_tests.cpp)
  target_link_libraries(dynsyslab_panel_workspace_tests PRIVATE dynsyslab_ui_support)
  add_test(NAME dynsyslab_panel_workspace_tests COMMAND dynsyslab_panel_workspace_tests)
  add_executable(dynsyslab_panel_catalog_tests tests/panel_catalog_tests.cpp)
  target_link_libraries(dynsyslab_panel_catalog_tests PRIVATE dynsyslab_viewer_services)
  add_test(NAME dynsyslab_panel_catalog_tests COMMAND dynsyslab_panel_catalog_tests)
  add_executable(dynsyslab_plot_interaction_tests tests/plot_interaction_tests.cpp)
  target_link_libraries(dynsyslab_plot_interaction_tests PRIVATE dynsyslab_viewer_interaction)
  add_test(NAME dynsyslab_plot_interaction_tests COMMAND dynsyslab_plot_interaction_tests)
  set_tests_properties(dynsyslab_plot_interaction_tests PROPERTIES
    LABELS "native;ui;plot-interaction")
  add_executable(dynsyslab_plot_3d_tests tests/plot_3d_tests.cpp)
  target_link_libraries(dynsyslab_plot_3d_tests PRIVATE dynsyslab_viewer_interaction)
  add_test(NAME dynsyslab_plot_3d_tests COMMAND dynsyslab_plot_3d_tests)
  set_tests_properties(dynsyslab_plot_3d_tests PROPERTIES
    LABELS "native;ui;plot-interaction;three-dimensional")
  add_executable(dynsyslab_plot_3d_canvas_tests tests/plot_3d_canvas_tests.cpp)
  target_link_libraries(dynsyslab_plot_3d_canvas_tests PRIVATE dynsyslab_ui_support)
  add_test(NAME dynsyslab_plot_3d_canvas_tests COMMAND dynsyslab_plot_3d_canvas_tests)
  set_tests_properties(dynsyslab_plot_3d_canvas_tests PROPERTIES
    LABELS "native;ui;plot-interaction;three-dimensional")
  add_executable(dynsyslab_plot_cursor_tests tests/plot_cursor_tests.cpp)
  target_link_libraries(dynsyslab_plot_cursor_tests PRIVATE dynsyslab_ui_support)
  add_test(NAME dynsyslab_plot_cursor_tests COMMAND dynsyslab_plot_cursor_tests)
  set_tests_properties(dynsyslab_plot_cursor_tests PROPERTIES
    LABELS "native;ui;plot-interaction")
  add_test(NAME dynsyslab_viewer_cli_interaction_contract
    COMMAND dynsyslab_viewer_cli verify-interactions --json)
  set_tests_properties(dynsyslab_viewer_cli_interaction_contract PROPERTIES
    PASS_REGULAR_EXPRESSION "\"passed\":true"
    LABELS "native;ui;plot-interaction;cli")
endif()
