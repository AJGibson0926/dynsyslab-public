file(MAKE_DIRECTORY "${TEST_ROOT}/archives" "${TEST_ROOT}/cache")
file(WRITE "${TEST_ROOT}/plot.json" [=[
{"schema":"dynsyslab.viewer.payload/v1","view":"phase","tile":{"name":"CLI fixture","columns":1},
 "figures":[{"kind":"phase","title":"Phase","x_label":"r","y_label":"dr/dtau","coordinates":"dimensionless",
 "color_label":"\\tau","square":false,"traces":[{"label":"solution","x":[1,2,1],"y":[0,1,0],
 "encoding":"time","color_values":[0,1,2]}],
 "references":[{"x":1,"y":0,"label":"equilibrium","equilibrium":true}],
 "alternate_coordinates":{"name":"dimensional","x_label":"R (m)","y_label":"dR/dt (m/s)",
 "x_scale":0.00001,"y_scale":3.33333333333333,"time_scale":0.000003}}]}
]=])
file(WRITE "${TEST_ROOT}/catalog.json" [=[
{"schema":"dynsyslab.viewer.catalog/v1","items":[
 {"id":"parent","label":"Parent","kind":"cluster","source":"parent.json","category":"Training data"},
 {"id":"visible","parent":"parent","label":"Selected trajectory","kind":"trajectory","source":"plot.json","views":["phase"]},
 {"id":"masked","label":"Masked trajectory","kind":"trajectory","source":"masked.json","views":["phase"]}]}
]=])
file(WRITE "${TEST_ROOT}/focus.json" [=[
{"schema":"dynsyslab.viewer.archive-focus/v1","active":true,"label":"CLI test Mask","allowed_sources":["plot.json"]}
]=])
set(common --project "${TEST_ROOT}" --archives "${TEST_ROOT}/archives" --cache "${TEST_ROOT}/cache"
  --catalog "${TEST_ROOT}/catalog.json" --focus "${TEST_ROOT}/focus.json")
execute_process(COMMAND "${VIEWER_CLI}" tree ${common} RESULT_VARIABLE status OUTPUT_VARIABLE tree ERROR_VARIABLE error)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "Background tree failed: ${error}")
endif()
string(JSON visible GET "${tree}" visible_item_count)
string(JSON hidden GET "${tree}" hidden_item_count)
if(NOT visible EQUAL 2 OR NOT hidden EQUAL 1)
  message(FATAL_ERROR "Mask or parent retention failed: ${tree}")
endif()
execute_process(COMMAND "${VIEWER_CLI}" tree ${common} --parent parent --search selected
  RESULT_VARIABLE status OUTPUT_VARIABLE children)
string(JSON count GET "${children}" returned_item_count)
if(NOT status EQUAL 0 OR NOT count EQUAL 1)
  message(FATAL_ERROR "Data-tree child navigation/search failed")
endif()
execute_process(COMMAND "${VIEWER_CLI}" inspect ${common} --id masked RESULT_VARIABLE status ERROR_QUIET)
if(status EQUAL 0)
  message(FATAL_ERROR "CLI selected a masked object")
endif()
execute_process(COMMAND "${VIEWER_CLI}" inspect ${common} --source plot.json --view phase
  --coordinates dimensional --samples RESULT_VARIABLE status OUTPUT_VARIABLE plot ERROR_VARIABLE error)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "Plot inspection failed: ${error}")
endif()
string(JSON radius GET "${plot}" figures 0 traces 0 x 0)
string(JSON coordinate GET "${plot}" figures 0 coordinates)
string(JSON reference GET "${plot}" figures 0 references 0 x)
if(NOT coordinate STREQUAL "dimensional" OR NOT radius EQUAL 0.00001 OR NOT reference EQUAL 0.00001)
  message(FATAL_ERROR "Coordinate conversion or reference metadata was lost")
endif()
# Round-trip the exported full payload through the same loader.
file(WRITE "${TEST_ROOT}/plot.json" "${plot}")
execute_process(COMMAND "${VIEWER_CLI}" inspect ${common} --id visible --view phase
  --coordinates dimensionless --samples RESULT_VARIABLE status OUTPUT_VARIABLE roundtrip)
string(JSON radius GET "${roundtrip}" figures 0 traces 0 x 0)
if(NOT status EQUAL 0 OR radius LESS 0.999999999 OR radius GREATER 1.000000001)
  message(FATAL_ERROR "Exported plot did not round trip")
endif()
execute_process(COMMAND "${VIEWER_CLI}" render ${common} --id visible --view phase
  --output "${TEST_ROOT}/phase.png" --width 900 --height 500 --figure 0
  RESULT_VARIABLE status OUTPUT_VARIABLE rendered ERROR_VARIABLE error)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "Background rendering failed: ${error}")
endif()
string(JSON windows GET "${rendered}" visible_windows)
file(READ "${TEST_ROOT}/phase.png" signature LIMIT 8 HEX)
if(NOT windows EQUAL 0 OR NOT signature STREQUAL "89504e470d0a1a0a")
  message(FATAL_ERROR "Render must produce a PNG without visible windows")
endif()
execute_process(COMMAND "${VIEWER_CLI}" render ${common} --id visible --view phase
  --output "${TEST_ROOT}/bad.png" --figure 99 RESULT_VARIABLE status ERROR_QUIET)
if(status EQUAL 0)
  message(FATAL_ERROR "CLI accepted an invalid figure")
endif()
message(STATUS "Background tree, Mask, selection, inspection, coordinates, payload round trip, PNG: PASS")
