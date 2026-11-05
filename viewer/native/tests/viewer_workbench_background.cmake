file(MAKE_DIRECTORY "${TEST_ROOT}/archives" "${TEST_ROOT}/cache")
file(WRITE "${TEST_ROOT}/plot.json" [=[
{"schema":"dynsyslab.viewer.payload/v1","view":"phase","tile":{"name":"Workbench fixture","columns":1},
 "figures":[{"kind":"phase","title":"Phase","x_label":"r","y_label":"dr/dtau",
 "traces":[{"label":"solution","x":[1,2,1],"y":[0,1,0]}]}]}
]=])
file(WRITE "${TEST_ROOT}/catalog.json" [=[
{"schema":"dynsyslab.viewer.catalog/v1","items":[
 {"id":"kept","label":"Kept trajectory","kind":"trajectory","source":"plot.json","views":["phase"]},
 {"id":"hidden","label":"Hidden trajectory","kind":"trajectory","source":"hidden.json","views":["phase"]}]}
]=])
set(common --background-check --project-root "${TEST_ROOT}" --archive-root "${TEST_ROOT}/archives"
  --catalog "${TEST_ROOT}/catalog.json" --cache "${TEST_ROOT}/cache" --archive-focus "${TEST_ROOT}/focus.json")
file(WRITE "${TEST_ROOT}/focus.json" [=[
{"schema":"dynsyslab.viewer.archive-focus/v1","active":true,"label":"Desktop Mask","allowed_sources":["plot.json"]}
]=])
execute_process(COMMAND "${VIEWER}" ${common} --background-refresh-count 3
  RESULT_VARIABLE status TIMEOUT 40)
file(READ "${TEST_ROOT}/cache/workbench-check.json" report)
string(JSON count GET "${report}" catalog_items)
string(JSON nodes GET "${report}" navigation_nodes)
string(JSON label GET "${report}" focus_label)
string(JSON visible GET "${report}" visible)
string(JSON plot GET "${report}" has_plot)
string(JSON refreshes GET "${report}" refreshes)
string(JSON sorting GET "${report}" tree_sort_controls_checked)
if(NOT status EQUAL 0 OR NOT count EQUAL 1 OR NOT nodes EQUAL 2 OR visible
   OR NOT plot OR NOT refreshes EQUAL 3 OR NOT sorting OR
   NOT label STREQUAL "Desktop Mask")
  message(FATAL_ERROR "Hidden workbench startup, painting and repeated refresh failed: ${report}")
endif()
file(WRITE "${TEST_ROOT}/focus.json" [=[
{"schema":"dynsyslab.viewer.archive-focus/v1","active":true,"label":"Desktop Mask","allowed_sources":["absent.json"]}
]=])
execute_process(COMMAND "${VIEWER}" ${common} RESULT_VARIABLE status TIMEOUT 40)
file(READ "${TEST_ROOT}/cache/workbench-check.json" report)
string(JSON count GET "${report}" catalog_items)
string(JSON label GET "${report}" focus_label)
if(NOT status EQUAL 0 OR NOT count EQUAL 0 OR NOT label STREQUAL "Desktop Mask")
  message(FATAL_ERROR "An unmatched desktop Mask must remain active: ${report}")
endif()
file(WRITE "${TEST_ROOT}/focus.json" "{")
execute_process(COMMAND "${VIEWER}" ${common} RESULT_VARIABLE status TIMEOUT 40)
file(READ "${TEST_ROOT}/cache/workbench-check.json" report)
string(JSON count GET "${report}" catalog_items)
string(JSON error GET "${report}" error)
if(status EQUAL 0 OR NOT count EQUAL 0 OR NOT error MATCHES "Viewer Mask")
  message(FATAL_ERROR "An unreadable desktop Mask must not expose the full tree: ${report}")
endif()
message(STATUS "Hidden desktop Mask, painting, repeated refresh, unmatched and invalid Masks: PASS")

file(MAKE_DIRECTORY "${TEST_ROOT}/cache-unmasked")
file(WRITE "${TEST_ROOT}/catalog-unmasked.json" [=[
{"schema":"dynsyslab.viewer.catalog/v1","items":[
 {"id":"other","label":"Other data","kind":"trajectory","category":"Data","system":"Fixture","source":"plot.json","views":["phase"],"modified":100},
 {"id":"training","label":"Training data","kind":"presentation_group","category":"Data","system":"Fixture"},
 {"id":"training-old","label":"Training set #1","kind":"trajectory","category":"Data","system":"Fixture","parent":"training","source":"plot.json","views":["phase"],"modified":10},
 {"id":"training-new","label":"Training set #2","kind":"trajectory","category":"Data","system":"Fixture","parent":"training","source":"plot.json","views":["phase"],"modified":20}]}
]=])
execute_process(COMMAND "${VIEWER}" --background-check --project-root "${TEST_ROOT}"
  --archive-root "${TEST_ROOT}/archives" --catalog "${TEST_ROOT}/catalog-unmasked.json"
  --cache "${TEST_ROOT}/cache-unmasked" RESULT_VARIABLE status TIMEOUT 40)
file(READ "${TEST_ROOT}/cache-unmasked/workbench-check.json" report)
string(JSON count GET "${report}" catalog_items)
string(JSON label GET "${report}" focus_label)
string(JSON hidden GET "${report}" hidden_items)
string(JSON focus_path GET "${report}" focus_path)
string(JSON selected GET "${report}" selected_id)
string(JSON training_selected GET "${report}" training_data_selected)
string(JSON plot GET "${report}" has_plot)
if(NOT status EQUAL 0 OR NOT count EQUAL 4 OR NOT label STREQUAL "" OR
   NOT hidden EQUAL 0 OR NOT focus_path STREQUAL "" OR
   NOT selected STREQUAL "training-new" OR NOT training_selected OR NOT plot)
  message(FATAL_ERROR "Default startup did not show the full catalog at training data: ${report}")
endif()
message(STATUS "Unmasked default startup and newest training-data selection: PASS")

execute_process(COMMAND "${OBSERVABLE_FIXTURE}" --write-fixture "${TEST_ROOT}/trajectory.npz"
  RESULT_VARIABLE status)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "Could not write observable trajectory fixture")
endif()
file(WRITE "${TEST_ROOT}/catalog.json" [=[
{"schema":"dynsyslab.viewer.catalog/v1","items":[
 {"id":"trajectory","label":"Main trajectory","kind":"trajectory","source":"trajectory.npz","views":["phase"]}]}
]=])
file(WRITE "${TEST_ROOT}/focus.json" [=[
{"schema":"dynsyslab.viewer.archive-focus/v1","active":true,"label":"Desktop Mask","allowed_sources":["trajectory.npz"]}
]=])
execute_process(COMMAND "${VIEWER}" ${common} --background-workspace-check
  RESULT_VARIABLE status TIMEOUT 40)
file(READ "${TEST_ROOT}/cache/workbench-check.json" report)
string(JSON tabs GET "${report}" workspace_tabs_checked)
string(JSON figures GET "${report}" observable_figures)
string(JSON visible GET "${report}" visible)
# Eight Duffing scalar outputs, including intrinsic force, each have two tiles.
if(NOT status EQUAL 0 OR NOT tabs OR NOT figures EQUAL 16 OR visible)
  message(FATAL_ERROR "Main/Observables tab controls or state retention failed: ${report}")
endif()
execute_process(COMMAND "${VIEWER_CLI}" inspect --id trajectory --workspace observables
  --project "${TEST_ROOT}" --archives "${TEST_ROOT}/archives" --cache "${TEST_ROOT}/cache"
  --catalog "${TEST_ROOT}/catalog.json" --focus "${TEST_ROOT}/focus.json"
  OUTPUT_VARIABLE inspection RESULT_VARIABLE status)
string(JSON figures LENGTH "${inspection}" figures)
if(NOT status EQUAL 0 OR NOT figures EQUAL 16)
  message(FATAL_ERROR "CLI observable workspace differs from the desktop: ${inspection}")
endif()
message(STATUS "Hidden workspace tabs, retained Main, complete observables and CLI parity: PASS")
