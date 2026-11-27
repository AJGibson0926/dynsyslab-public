get_filename_component(DYNSYSLAB_ROOT "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)
set(DYNSYSLAB_VIEWER_NATIVE "${DYNSYSLAB_ROOT}/viewer/native")

add_library(dynsyslab_viewer_interaction STATIC
  ${DYNSYSLAB_VIEWER_NATIVE}/src/MathText.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/PlotScale.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/PlotPresentation.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/Plot3D.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/PlotInteractionModel.cpp)
target_include_directories(dynsyslab_viewer_interaction PUBLIC
  ${DYNSYSLAB_VIEWER_NATIVE}/include
  ${DYNSYSLAB_ROOT}/native/include)
target_compile_features(dynsyslab_viewer_interaction PUBLIC cxx_std_17)
if(MSVC)
  target_compile_options(dynsyslab_viewer_interaction PRIVATE
    /W4 /permissive- /EHsc /utf-8)
else()
  target_compile_options(dynsyslab_viewer_interaction PRIVATE
    -Wall -Wextra -Wpedantic)
endif()

add_library(dynsyslab_ui_support STATIC
  ${DYNSYSLAB_VIEWER_NATIVE}/src/Model.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/PanelWorkspace.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/CatalogPersistence.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/CatalogDiscovery.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/CatalogPresentation.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/CatalogPresentationText.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/CatalogTreeSort.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/GdiBackBuffer.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/PlotFigureDrawing.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/Plot3DDrawing.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/PlotCanvasScrolling.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/PlotCanvasRaster.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/PlotCanvas.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/PlotCanvasCoordinates.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/PlotCanvasState.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/FolderTab.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/PlotSeries.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/PlotHeatmap.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/PlotColorbar.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/ComplexHeatmap.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/MathTextDrawing.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/PlotCanvasInteraction.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/PlotCanvas3DInteraction.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/PlotCanvasComposition.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/PlotKromSelection.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/PlotDragToolbar.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/PlotLegend.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/PlotProbeOverlay.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/PlotPanelChrome.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/ArchiveCatalog.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/SystemViewDefaults.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/ArchiveFocus.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/HierarchyCatalog.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/Process.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/RuntimePaths.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/StartupDiagnostics.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/Theme.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/ToggleButton.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/HorizontalPaneLayout.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/VerticalSplitter.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/OrnamentalScrollBar.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/OrnamentalComboBox.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/NavigationPane.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/InspectorPane.cpp
  ${DYNSYSLAB_VIEWER_NATIVE}/src/Help.cpp)

target_include_directories(dynsyslab_ui_support PUBLIC
  ${DYNSYSLAB_VIEWER_NATIVE}/include
  ${DYNSYSLAB_ROOT}/native_ui/windows)
target_compile_features(dynsyslab_ui_support PUBLIC cxx_std_17)
target_compile_definitions(dynsyslab_ui_support PUBLIC UNICODE _UNICODE
  WIN32_LEAN_AND_MEAN NOMINMAX)
target_link_libraries(dynsyslab_ui_support PUBLIC comctl32 gdiplus shell32
  shlwapi dwmapi uxtheme dynsyslab_viewer_interaction
  dynsyslab_foundation dynsyslab_operations)
if(MSVC)
  target_compile_options(dynsyslab_ui_support PRIVATE /W4 /permissive- /EHsc /utf-8)
else()
  target_compile_options(dynsyslab_ui_support PRIVATE -Wall -Wextra -Wpedantic)
endif()

install(TARGETS dynsyslab_viewer_interaction ARCHIVE DESTINATION lib)
install(FILES
  ${DYNSYSLAB_VIEWER_NATIVE}/include/dynsyslab_viewer/Json.hpp
  ${DYNSYSLAB_VIEWER_NATIVE}/include/dynsyslab_viewer/Model.hpp
  ${DYNSYSLAB_VIEWER_NATIVE}/include/dynsyslab_viewer/MathText.hpp
  ${DYNSYSLAB_VIEWER_NATIVE}/include/dynsyslab_viewer/PlotScale.hpp
  ${DYNSYSLAB_VIEWER_NATIVE}/include/dynsyslab_viewer/PlotPresentation.hpp
  ${DYNSYSLAB_VIEWER_NATIVE}/include/dynsyslab_viewer/PlotInteractionModel.hpp
  ${DYNSYSLAB_VIEWER_NATIVE}/include/dynsyslab_viewer/Plot3D.hpp
  ${DYNSYSLAB_VIEWER_NATIVE}/include/dynsyslab_viewer/CatalogTreeSort.hpp
  DESTINATION include/dynsyslab_viewer)
