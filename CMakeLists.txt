cmake_minimum_required(VERSION 3.24)
project(DynSysLab VERSION 0.5.1 LANGUAGES CXX)

include(CTest)
include(cmake/DynSysLabExtensions.cmake)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
if(MSVC)
  set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
endif()

option(DYNSYSLAB_BUILD_TESTS "Build native release gates" ON)
option(DYNSYSLAB_BUILD_PYTHON "Build optional Python bindings" OFF)
option(DYNSYSLAB_BUILD_VIEWER "Build DynSysLab Viewer" ON)
option(DYNSYSLAB_BUILD_COMPUTE "Build DynSysLab Solver and Post-processor" ON)

set(DYNSYSLAB_BUILD_CLI OFF CACHE BOOL "" FORCE)
set(DYNSYSLAB_ENABLE_ARCHIVE_IO ON CACHE BOOL "" FORCE)
set(DYNSYSLAB_ENABLE_LINALG ON CACHE BOOL "" FORCE)
set(DYNSYSLAB_ENABLE_SUNDIALS ON CACHE BOOL "" FORCE)
set(DYNSYSLAB_ENABLE_SQLITE ON CACHE BOOL "" FORCE)
set(DYNSYSLAB_BUILD_SHARED_MODULES ON CACHE BOOL "" FORCE)
set(DYNSYSLAB_COMPUTE_BUILD_TESTS ${DYNSYSLAB_BUILD_TESTS} CACHE BOOL "" FORCE)
set(DYNSYSLAB_VIEWER_BUILD_TESTS ${DYNSYSLAB_BUILD_TESTS} CACHE BOOL "" FORCE)

add_subdirectory(native)
add_subdirectory(workflow/native)
add_subdirectory(native_ui)
if(DYNSYSLAB_BUILD_COMPUTE)
  add_subdirectory(compute_workbench/native)
endif()
if(DYNSYSLAB_BUILD_VIEWER)
  add_subdirectory(viewer/native)
endif()

if(DYNSYSLAB_BUILD_TESTS)
  # Discover test executables after all subdirectories are registered. A hand-maintained
  # list can silently run stale ABI-incompatible binaries after public header changes.
  function(_dynsyslab_collect_test_targets directory output)
    get_property(_targets DIRECTORY "${directory}" PROPERTY BUILDSYSTEM_TARGETS)
    set(_tests)
    foreach(_target IN LISTS _targets)
      if(_target MATCHES "^dynsyslab_.*(_tests|_test_extension)$")
        list(APPEND _tests "${_target}")
      endif()
    endforeach()
    get_property(_directories DIRECTORY "${directory}" PROPERTY SUBDIRECTORIES)
    foreach(_directory IN LISTS _directories)
      _dynsyslab_collect_test_targets("${_directory}" _child_tests)
      list(APPEND _tests ${_child_tests})
    endforeach()
    set(${output} ${_tests} PARENT_SCOPE)
  endfunction()
  _dynsyslab_collect_test_targets("${CMAKE_CURRENT_SOURCE_DIR}" _dsl_test_targets)
  # Scripted CTest contracts invoke the packaged applications and CLIs as well.
  foreach(_application dynsyslab_compute dynsyslab_viewer dynsyslab_cli dynsyslab_viewer_cli)
    if(TARGET ${_application})
      list(APPEND _dsl_test_targets ${_application})
    endif()
  endforeach()
  add_custom_target(dynsyslab_test_binaries)
  add_dependencies(dynsyslab_test_binaries ${_dsl_test_targets})
  unset(_dsl_test_targets)
endif()

if(WIN32)
  set(CMAKE_INSTALL_SYSTEM_RUNTIME_DESTINATION bin)
  set(CMAKE_INSTALL_UCRT_LIBRARIES TRUE)
  include(InstallRequiredSystemLibraries)
endif()

install(DIRECTORY docs/ DESTINATION resources/docs
  PATTERN "archive" EXCLUDE
  PATTERN "assets" EXCLUDE
  PATTERN "DISK_SPACE_AUDIT*" EXCLUDE
  PATTERN "*.pdf" EXCLUDE
  PATTERN "*.png" EXCLUDE
  PATTERN "*.jpg" EXCLUDE
  PATTERN "__pycache__" EXCLUDE
  PATTERN "*.pyc" EXCLUDE
  PATTERN "*.py" EXCLUDE)
install(FILES config/antenna_policy.template.json
              config/antenna_suite.template.json
              config/s3_snapshot.example.json
  DESTINATION resources/config)
install(FILES README.md ARCHITECTURE.md SOLVER.md PERFORMANCE.md CHANGELOG.md
  DESTINATION resources)
install(FILES install.ps1 DESTINATION .)
install(DIRECTORY examples/native/ DESTINATION resources/examples
  PATTERN "case-studies" EXCLUDE)
if(EXISTS "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/share")
  install(DIRECTORY "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/share/"
    DESTINATION resources/third-party-licenses FILES_MATCHING PATTERN "copyright")
endif()
install(FILES native_ui/assets/dynsyslab-app.ico native_ui/assets/dynsyslab-viewer.ico DESTINATION resources/assets)
install(DIRECTORY examples/native-extension/
  DESTINATION resources/sdk/examples/native-extension)

if(TARGET dynsyslab_compute)
  install(FILES $<TARGET_RUNTIME_DLLS:dynsyslab_compute> DESTINATION bin)
endif()
if(TARGET dynsyslab_viewer)
  install(FILES $<TARGET_RUNTIME_DLLS:dynsyslab_viewer> DESTINATION bin)
endif()
if(WIN32 AND TARGET dynsyslab_cli)
  set(_dsl_runtime_search "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/bin")
  install(CODE "
    file(GET_RUNTIME_DEPENDENCIES
      EXECUTABLES \"$<TARGET_FILE:dynsyslab_cli>\"
      RESOLVED_DEPENDENCIES_VAR _dsl_dependencies
      UNRESOLVED_DEPENDENCIES_VAR _dsl_unresolved
      DIRECTORIES \"${_dsl_runtime_search}\"
      PRE_EXCLUDE_REGEXES \"api-ms-.*\" \"ext-ms-.*\")
    file(TO_CMAKE_PATH \"${_dsl_runtime_search}/\" _dsl_runtime_root)
    string(TOLOWER \"\${_dsl_runtime_root}\" _dsl_runtime_root)
    file(TO_CMAKE_PATH \"${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/\" _dsl_build_runtime_root)
    string(TOLOWER \"\${_dsl_build_runtime_root}\" _dsl_build_runtime_root)
    foreach(_dsl_dependency IN LISTS _dsl_dependencies)
      file(TO_CMAKE_PATH \"\${_dsl_dependency}\" _dsl_dependency_path)
      string(TOLOWER \"\${_dsl_dependency_path}\" _dsl_dependency_path_lower)
      string(FIND \"\${_dsl_dependency_path_lower}\" \"\${_dsl_runtime_root}\" _dsl_runtime_index)
      string(FIND \"\${_dsl_dependency_path_lower}\" \"\${_dsl_build_runtime_root}\" _dsl_build_runtime_index)
      if(_dsl_runtime_index EQUAL 0 OR _dsl_build_runtime_index EQUAL 0)
        file(INSTALL DESTINATION \"\${CMAKE_INSTALL_PREFIX}/bin\"
          TYPE SHARED_LIBRARY FILES \"\${_dsl_dependency}\")
      endif()
    endforeach()
    foreach(_dsl_missing IN LISTS _dsl_unresolved)
      string(TOLOWER \"\${_dsl_missing}\" _dsl_missing_lower)
      if(_dsl_missing_lower MATCHES \"^(dynsyslab|sundials|sqlite|archive|liblzma|z).*\\\\.dll$\")
        message(FATAL_ERROR \"Unresolved required DynSysLab runtime dependency: \${_dsl_missing}\")
      endif()
    endforeach()
  ")
  unset(_dsl_runtime_search)
endif()

set(CPACK_PACKAGE_NAME "DynSysLab")
set(CPACK_PACKAGE_VENDOR "DynSysLab")
set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
set(CPACK_PACKAGE_FILE_NAME "DynSysLab-${PROJECT_VERSION}-windows-x64")
set(CPACK_GENERATOR "ZIP")
include(CPack)

if(DYNSYSLAB_BUILD_TESTS AND TARGET dynsyslab_compute AND TARGET dynsyslab_viewer)
  add_custom_target(dynsyslab_package_smoke
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}/package-smoke
    COMMAND ${CMAKE_COMMAND} --install ${CMAKE_BINARY_DIR} --config $<CONFIG>
      --prefix ${CMAKE_BINARY_DIR}/package-smoke
    COMMAND ${CMAKE_COMMAND}
      -DSTAGE_DIR=${CMAKE_BINARY_DIR}/package-smoke
      -P ${CMAKE_SOURCE_DIR}/cmake/PackageSmoke.cmake
    DEPENDS dynsyslab_compute dynsyslab_viewer dynsyslab_viewer_cli dynsyslab_cli dynsyslab_linalg
    USES_TERMINAL)
endif()
