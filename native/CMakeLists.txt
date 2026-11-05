cmake_minimum_required(VERSION 3.20)
project(dynsyslab_native VERSION 0.5.0 LANGUAGES CXX)

include(CTest)
option(DYNSYSLAB_BUILD_PYTHON "Build optional Python bindings" OFF)
option(DYNSYSLAB_BUILD_CLI "Build the small core diagnostic CLI" ON)
option(DYNSYSLAB_BUILD_TESTS "Build native-core tests" ON)
option(DYNSYSLAB_BUILD_SHARED_MODULES "Build versioned native DLL modules" ON)
option(DYNSYSLAB_ENABLE_ARCHIVE_IO "Build NPZ/LZMA archive compatibility" OFF)
option(DYNSYSLAB_ENABLE_LINALG "Build Eigen-backed scientific methods" OFF)
option(DYNSYSLAB_ENABLE_SUNDIALS "Build the SUNDIALS BDF backend" OFF)
option(DYNSYSLAB_ENABLE_SQLITE "Build the durable native job store" OFF)
set(DYNSYSLAB_EIGEN_COMPILE_JOBS 3 CACHE STRING
  "Maximum concurrent Eigen-heavy compiler processes")

if(MSVC)
  # STL values cross module boundaries, so every DLL shares one CRT heap.
  # Packaging installs the redistributable runtime beside the executables.
  set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
endif()
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
endif()

set(DYNSYSLAB_LIBRARY_KIND STATIC)
if(DYNSYSLAB_BUILD_SHARED_MODULES)
  set(DYNSYSLAB_LIBRARY_KIND SHARED)
endif()

function(dynsyslab_module target)
  target_include_directories(${target} PUBLIC include)
  target_compile_features(${target} PUBLIC cxx_std_17)
  set_target_properties(${target} PROPERTIES POSITION_INDEPENDENT_CODE ON)
  if(DYNSYSLAB_BUILD_SHARED_MODULES AND WIN32)
    set_target_properties(${target} PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
  endif()
  if(WIN32)
    target_compile_definitions(${target} PUBLIC WIN32_LEAN_AND_MEAN NOMINMAX)
  endif()
  if(MSVC)
    target_compile_options(${target} PRIVATE /W4 /permissive- /EHsc /utf-8 /MP4 /bigobj)
  else()
    target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
  endif()
endfunction()

function(dynsyslab_eigen_module target)
  dynsyslab_module(${target})
  target_link_libraries(${target} PRIVATE Eigen3::Eigen)
  target_precompile_headers(${target} PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/src/EigenPch.hpp")
  set_property(TARGET ${target} PROPERTY JOB_POOL_COMPILE dynsyslab_eigen_compile)
endfunction()

add_library(dynsyslab_foundation ${DYNSYSLAB_LIBRARY_KIND}
  src/Models.cpp src/Npy.cpp src/Json.cpp src/Text.cpp src/FileDigest.cpp
  src/LegacyCsv.cpp src/LegacyCsvImport.cpp
  src/TimeSeries.cpp src/Fourier.cpp src/DiscreteFourier.cpp src/TimeSeriesArchive.cpp
  src/Extension.cpp src/ODE.cpp src/Integrator.cpp src/Hierarchy.cpp
  src/AnalyticExpression.cpp src/ObservableExpressions.cpp
  src/SpecialFunctions.cpp src/Lauricella.cpp src/QuarticOrbit.cpp src/AnalyticSolution.cpp
  src/AnalyticQuadrature.cpp src/AnalyticDoublePendulum.cpp src/AnalyticModelLimits.cpp src/AnalyticKuramoto.cpp
  src/NeuralKnowledge.cpp src/CoupledKnowledge.cpp src/BubbleKnowledge.cpp
  src/AnalyticHenonHeiles.cpp src/AnalyticOscillators.cpp src/AnalyticRPE.cpp src/AnalyticOtherSystems.cpp
  src/SystemsHarmonic.cpp src/ResearchPapers.cpp src/ResearchPapersExtra.cpp src/DynamicalKnowledge.cpp src/RpeDynamicalKnowledge.cpp src/DoublePendulumKnowledge.cpp
  src/Observable.cpp src/System.cpp src/SystemSupport.cpp src/SystemLabels.cpp src/Targets.cpp
  src/TrainingDomain.cpp
  src/SystemsBasic.cpp src/SystemsMechanical.cpp src/SystemsBiological.cpp
  src/SystemsBubble.cpp src/SystemsRPE.cpp src/SystemsMarmottant.cpp
  src/SystemsReddySzeri.cpp)
dynsyslab_module(dynsyslab_foundation)
find_package(boost_math CONFIG REQUIRED)
find_package(boost_multiprecision CONFIG REQUIRED)
target_link_libraries(dynsyslab_foundation PRIVATE Boost::math)

add_library(dynsyslab_delay ${DYNSYSLAB_LIBRARY_KIND} src/Delay.cpp
  src/AnalyticDelay.cpp src/DelayExponential.cpp src/DelayResearch.cpp)
dynsyslab_module(dynsyslab_delay)
target_link_libraries(dynsyslab_delay PUBLIC dynsyslab_foundation PRIVATE Boost::math Boost::multiprecision)

add_library(dynsyslab_lattice ${DYNSYSLAB_LIBRARY_KIND} src/Lattice.cpp src/LatticeResearch.cpp)
dynsyslab_module(dynsyslab_lattice)
target_link_libraries(dynsyslab_lattice PUBLIC dynsyslab_foundation)

add_library(dynsyslab_pde ${DYNSYSLAB_LIBRARY_KIND} src/PDE.cpp src/PDEParallel.cpp src/AnalyticHeat.cpp src/HeatResearch.cpp)
dynsyslab_module(dynsyslab_pde)
target_link_libraries(dynsyslab_pde PUBLIC dynsyslab_foundation)

if(DYNSYSLAB_ENABLE_ARCHIVE_IO)
  find_package(LibArchive REQUIRED)
  if(MSVC AND DEFINED VCPKG_INSTALLED_DIR AND DEFINED VCPKG_TARGET_TRIPLET)
    set(_dsl_archive_root "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}")
    set_property(TARGET LibArchive::LibArchive PROPERTY IMPORTED_CONFIGURATIONS
      "DEBUG;RELEASE;RELWITHDEBINFO;MINSIZEREL")
    set_target_properties(LibArchive::LibArchive PROPERTIES
      IMPORTED_LOCATION_DEBUG "${_dsl_archive_root}/debug/lib/archive.lib"
      IMPORTED_LOCATION_RELEASE "${_dsl_archive_root}/lib/archive.lib"
      IMPORTED_LOCATION_RELWITHDEBINFO "${_dsl_archive_root}/lib/archive.lib"
      IMPORTED_LOCATION_MINSIZEREL "${_dsl_archive_root}/lib/archive.lib")
    unset(_dsl_archive_root)
  endif()
  target_sources(dynsyslab_foundation PRIVATE src/Npz.cpp)
  target_link_libraries(dynsyslab_foundation PUBLIC LibArchive::LibArchive)
  target_compile_definitions(dynsyslab_foundation PUBLIC DYNSYSLAB_ARCHIVE_IO=1)
else()
  target_sources(dynsyslab_foundation PRIVATE src/NpzUnavailable.cpp)
endif()

if(DYNSYSLAB_ENABLE_SUNDIALS)
  find_package(SUNDIALS CONFIG REQUIRED)
  target_sources(dynsyslab_foundation PRIVATE src/IntegratorSundials.cpp)
  if(TARGET SUNDIALS::cvode_shared)
    target_link_libraries(dynsyslab_foundation PUBLIC SUNDIALS::cvode_shared
      SUNDIALS::nvecserial_shared SUNDIALS::sunmatrixdense_shared
      SUNDIALS::sunlinsoldense_shared)
  else()
    target_link_libraries(dynsyslab_foundation PUBLIC SUNDIALS::cvode_static
      SUNDIALS::nvecserial_static SUNDIALS::sunmatrixdense_static
      SUNDIALS::sunlinsoldense_static)
  endif()
  target_compile_definitions(dynsyslab_foundation PUBLIC DYNSYSLAB_SUNDIALS=1)
else()
  target_sources(dynsyslab_foundation PRIVATE src/IntegratorSundialsUnavailable.cpp)
endif()

add_library(dynsyslab_eigenfunctions ${DYNSYSLAB_LIBRARY_KIND}
  src/EigenfunctionSpectrum.cpp src/EigenfunctionCoordinates.cpp
  src/PolynomialEigenfunctions.cpp src/EigenfunctionFit.cpp src/PolynomialBasis.cpp
  src/PolynomialEigenfunctionConstruction.cpp src/PolynomialEigenfunctionSpectrum.cpp)
dynsyslab_module(dynsyslab_eigenfunctions)
target_link_libraries(dynsyslab_eigenfunctions PUBLIC dynsyslab_foundation)

if(DYNSYSLAB_ENABLE_LINALG)
  find_package(Eigen3 CONFIG REQUIRED)
  if(NOT DYNSYSLAB_EIGEN_COMPILE_JOBS MATCHES "^[1-9][0-9]*$")
    message(FATAL_ERROR "DYNSYSLAB_EIGEN_COMPILE_JOBS must be a positive integer")
  endif()
  set_property(GLOBAL APPEND PROPERTY JOB_POOLS
    dynsyslab_eigen_compile=${DYNSYSLAB_EIGEN_COMPILE_JOBS})

  add_library(dynsyslab_numerics ${DYNSYSLAB_LIBRARY_KIND}
    src/EigenNumericsSvd.cpp src/EigenNumericsLeastSquaresReal.cpp
    src/EigenNumericsLeastSquaresComplex.cpp
    src/EigenNumericsNormalEquations.cpp
    src/EigenNumericsNullspace.cpp src/EigenNumericsComplexSpectrum.cpp
    src/EigenNumericsEigen.cpp src/EigenNumericsControl.cpp
    src/EigenNumericsControlFinite.cpp
    src/EigenNumericsMatrixFunctions.cpp
    src/EigenNumericsProductReal.cpp src/EigenNumericsProductComplex.cpp)
  dynsyslab_eigen_module(dynsyslab_numerics)

  add_library(dynsyslab_phase_space ${DYNSYSLAB_LIBRARY_KIND}
    src/PhaseSpace.cpp src/PhaseKnowledge.cpp src/PhaseFlow.cpp src/PhaseFixedPoints.cpp
    src/PhaseManifolds.cpp src/PhasePeriodicOrbits.cpp src/PhasePoincare.cpp
    src/PhaseBifurcation.cpp src/PhaseBifurcationEvents.cpp src/PhaseForcingClock.cpp
    src/PhaseForcedOrbits.cpp src/PhaseForcedManifolds.cpp
    src/PhaseSurfaces.cpp src/PhaseSurfaceGeometry.cpp src/PhaseSurfaceShooting.cpp
    src/PhaseSurfaceGrowth.cpp src/PhaseSurfaceValidation.cpp)
  dynsyslab_module(dynsyslab_phase_space)
  target_link_libraries(dynsyslab_phase_space PUBLIC dynsyslab_foundation PRIVATE dynsyslab_numerics)

  add_library(dynsyslab_identification ${DYNSYSLAB_LIBRARY_KIND}
    src/SINDy.cpp src/SindyEigenfunctionFit.cpp
    src/BoundaryEigenfunctions.cpp src/EigenvalueOptimization.cpp)
  dynsyslab_module(dynsyslab_identification)
  target_link_libraries(dynsyslab_identification PUBLIC dynsyslab_eigenfunctions
    PRIVATE dynsyslab_numerics)
else()
  add_library(dynsyslab_identification ${DYNSYSLAB_LIBRARY_KIND}
    src/SINDyUnavailable.cpp src/SindyEigenfunctionFit.cpp)
  dynsyslab_module(dynsyslab_identification)
  target_link_libraries(dynsyslab_identification PUBLIC dynsyslab_eigenfunctions)
endif()

add_library(dynsyslab_koopman ${DYNSYSLAB_LIBRARY_KIND}
  src/KROM.cpp src/KoopmanModel.cpp src/KoopmanModelConstruction.cpp src/KoopmanValidation.cpp
  src/ScatteredEigenfunctions.cpp src/KoopmanDatasets.cpp src/KoopmanInspection.cpp
  src/ScientificCompatibility.cpp)
dynsyslab_module(dynsyslab_koopman)
target_link_libraries(dynsyslab_koopman PUBLIC dynsyslab_foundation dynsyslab_eigenfunctions)
if(NOT DYNSYSLAB_ENABLE_LINALG)
  target_sources(dynsyslab_koopman PRIVATE src/SpectralUnavailable.cpp)
endif()
if(DYNSYSLAB_ENABLE_ARCHIVE_IO)
  target_sources(dynsyslab_koopman PRIVATE src/KoopmanArchive.cpp
    src/ControlMapArchive.cpp src/KoopmanPolynomialContract.cpp src/KoopmanBinding.cpp)
endif()

if(DYNSYSLAB_ENABLE_LINALG)
  add_library(dynsyslab_spectral ${DYNSYSLAB_LIBRARY_KIND}
    src/Spectral.cpp src/DmdControl.cpp src/DmdControlEigenfunctions.cpp src/EdmdControl.cpp)
  dynsyslab_module(dynsyslab_spectral)
  target_link_libraries(dynsyslab_spectral PUBLIC
    dynsyslab_foundation dynsyslab_koopman PRIVATE dynsyslab_numerics)
  if(DYNSYSLAB_ENABLE_ARCHIVE_IO)
    target_sources(dynsyslab_spectral PRIVATE src/EdmdControlArchive.cpp)
  endif()

  add_library(dynsyslab_decomposition ${DYNSYSLAB_LIBRARY_KIND}
    src/Decomposition.cpp)
  dynsyslab_module(dynsyslab_decomposition)
  target_link_libraries(dynsyslab_decomposition PUBLIC dynsyslab_foundation
    PRIVATE dynsyslab_numerics)

  add_library(dynsyslab_observation ${DYNSYSLAB_LIBRARY_KIND}
    src/HistoryProjection.cpp src/FrozenHistoryTransform.cpp)
  dynsyslab_eigen_module(dynsyslab_observation)
  target_link_libraries(dynsyslab_observation PUBLIC dynsyslab_foundation)
  if(DYNSYSLAB_ENABLE_ARCHIVE_IO)
    target_sources(dynsyslab_observation PRIVATE src/HistoryObservationArchive.cpp)
  endif()

  add_library(dynsyslab_history_koopman ${DYNSYSLAB_LIBRARY_KIND}
    src/HistoryKoopman.cpp src/HistoryInputRegression.cpp)
  dynsyslab_eigen_module(dynsyslab_history_koopman)
  target_link_libraries(dynsyslab_history_koopman PUBLIC
    dynsyslab_koopman dynsyslab_observation PRIVATE dynsyslab_numerics)

  add_library(dynsyslab_reconstruction ${DYNSYSLAB_LIBRARY_KIND}
    src/Reconstruction.cpp)
  dynsyslab_module(dynsyslab_reconstruction)
  target_link_libraries(dynsyslab_reconstruction PUBLIC dynsyslab_koopman
    PRIVATE dynsyslab_numerics)

  add_library(dynsyslab_pca ${DYNSYSLAB_LIBRARY_KIND} src/PCA.cpp)
  dynsyslab_module(dynsyslab_pca)
  target_link_libraries(dynsyslab_pca PUBLIC dynsyslab_spectral)

  add_library(dynsyslab_linalg SHARED src/LinalgCApi.cpp)
  target_include_directories(dynsyslab_linalg PUBLIC include)
  target_compile_features(dynsyslab_linalg PUBLIC cxx_std_17)
  target_compile_definitions(dynsyslab_linalg PRIVATE DYNSYSLAB_LINALG_BUILD=1)
  target_link_libraries(dynsyslab_linalg PRIVATE dynsyslab_numerics)
  if(MSVC)
    target_compile_options(dynsyslab_linalg PRIVATE /W4 /permissive- /EHsc /utf-8)
  endif()
  add_library(dynsyslab_control_maps ${DYNSYSLAB_LIBRARY_KIND}
    src/ControlRegression.cpp)
  dynsyslab_module(dynsyslab_control_maps)
  target_link_libraries(dynsyslab_control_maps PUBLIC dynsyslab_koopman
    PRIVATE dynsyslab_numerics)

  add_library(dynsyslab_control_screening ${DYNSYSLAB_LIBRARY_KIND}
    src/KromScreening.cpp)
  dynsyslab_module(dynsyslab_control_screening)
  target_link_libraries(dynsyslab_control_screening PUBLIC
    dynsyslab_koopman dynsyslab_spectral)

  add_library(dynsyslab_estimation ${DYNSYSLAB_LIBRARY_KIND}
    src/Estimation.cpp)
  dynsyslab_module(dynsyslab_estimation)
  target_link_libraries(dynsyslab_estimation PUBLIC
    dynsyslab_koopman PRIVATE dynsyslab_numerics)

  add_library(dynsyslab_control ${DYNSYSLAB_LIBRARY_KIND}
    src/LQR.cpp src/ClassicalControl.cpp)
  dynsyslab_module(dynsyslab_control)
  target_link_libraries(dynsyslab_control PUBLIC
    dynsyslab_koopman PRIVATE dynsyslab_numerics)
endif()

if(DYNSYSLAB_ENABLE_SQLITE)
  find_package(unofficial-sqlite3 CONFIG REQUIRED)
  add_library(dynsyslab_operations ${DYNSYSLAB_LIBRARY_KIND}
    src/JobStore.cpp src/ArchiveAudit.cpp src/ArchiveDeduplication.cpp
    src/ArchiveIndex.cpp src/ArchiveCatalogSnapshot.cpp)
  dynsyslab_module(dynsyslab_operations)
  target_link_libraries(dynsyslab_operations PUBLIC dynsyslab_foundation
    unofficial::sqlite3::sqlite3)
  target_compile_definitions(dynsyslab_operations PUBLIC DYNSYSLAB_SQLITE=1)
endif()

add_library(dynsyslab_core INTERFACE)
target_link_libraries(dynsyslab_core INTERFACE dynsyslab_foundation dynsyslab_delay
  dynsyslab_lattice dynsyslab_pde dynsyslab_koopman dynsyslab_identification)
if(TARGET dynsyslab_control_maps)
  target_link_libraries(dynsyslab_core INTERFACE dynsyslab_control_maps dynsyslab_control_screening)
endif()
if(TARGET dynsyslab_phase_space)
  target_link_libraries(dynsyslab_core INTERFACE dynsyslab_phase_space)
endif()
if(TARGET dynsyslab_eigenfunctions)
  target_link_libraries(dynsyslab_core INTERFACE dynsyslab_eigenfunctions)
endif()
if(TARGET dynsyslab_spectral)
  target_link_libraries(dynsyslab_core INTERFACE dynsyslab_spectral)
endif()
if(TARGET dynsyslab_decomposition)
  target_link_libraries(dynsyslab_core INTERFACE dynsyslab_decomposition)
endif()
if(TARGET dynsyslab_history_koopman)
  target_link_libraries(dynsyslab_core INTERFACE dynsyslab_history_koopman)
endif()
if(TARGET dynsyslab_reconstruction)
  target_link_libraries(dynsyslab_core INTERFACE dynsyslab_reconstruction)
endif()
if(TARGET dynsyslab_pca)
  target_link_libraries(dynsyslab_core INTERFACE dynsyslab_pca)
endif()
if(TARGET dynsyslab_control)
  target_link_libraries(dynsyslab_core INTERFACE dynsyslab_control)
endif()
if(TARGET dynsyslab_estimation)
  target_link_libraries(dynsyslab_core INTERFACE dynsyslab_estimation)
endif()

add_library(dynsyslab_c_api SHARED src/CApi.cpp)
target_include_directories(dynsyslab_c_api PUBLIC include)
target_compile_features(dynsyslab_c_api PUBLIC cxx_std_17)
target_compile_definitions(dynsyslab_c_api PRIVATE DYNSYSLAB_C_API_BUILD=1)
target_link_libraries(dynsyslab_c_api PRIVATE dynsyslab_foundation dynsyslab_delay
  dynsyslab_lattice dynsyslab_pde)
if(MSVC)
  target_compile_options(dynsyslab_c_api PRIVATE /W4 /permissive- /EHsc /utf-8)
endif()
if(DYNSYSLAB_ENABLE_LINALG AND DYNSYSLAB_ENABLE_ARCHIVE_IO)
    add_executable(dynsyslab_echo_input_experiment
      tools/echo_input/Main.cpp tools/echo_input/Experiment.cpp tools/echo_input/Refine.cpp
      tools/echo_input/EigenfunctionSweep.cpp)
    target_link_libraries(dynsyslab_echo_input_experiment PRIVATE dynsyslab_core)
endif()
if(DYNSYSLAB_BUILD_CLI)
  add_executable(dynsyslab_native_cli tools/main.cpp)
  target_link_libraries(dynsyslab_native_cli PRIVATE dynsyslab_core)
endif()

if(BUILD_TESTING AND DYNSYSLAB_BUILD_TESTS)
  foreach(test core system integrator analytic_solution)
    add_executable(dynsyslab_${test}_tests tests/${test}_tests.cpp)
    target_link_libraries(dynsyslab_${test}_tests PRIVATE dynsyslab_core)
    add_test(NAME dynsyslab_${test}_tests COMMAND dynsyslab_${test}_tests)
  endforeach()
  add_executable(dynsyslab_time_series_tests tests/time_series_tests.cpp)
  target_link_libraries(dynsyslab_time_series_tests PRIVATE dynsyslab_foundation)
  add_test(NAME dynsyslab_time_series_tests COMMAND dynsyslab_time_series_tests)
  set_tests_properties(dynsyslab_time_series_tests PROPERTIES
    LABELS "native;time-series;fourier;archive")
  add_executable(dynsyslab_force_observable_tests tests/force_observable_tests.cpp)
  target_link_libraries(dynsyslab_force_observable_tests PRIVATE dynsyslab_foundation)
  add_test(NAME dynsyslab_force_observable_tests COMMAND dynsyslab_force_observable_tests)
  add_executable(dynsyslab_catalogue_limits_tests tests/catalogue_limits_tests.cpp)
  target_link_libraries(dynsyslab_catalogue_limits_tests PRIVATE dynsyslab_foundation)
  add_test(NAME dynsyslab_catalogue_limits_tests COMMAND dynsyslab_catalogue_limits_tests)
  add_executable(dynsyslab_double_pendulum_analytic_tests tests/double_pendulum_analytic_tests.cpp)
  target_link_libraries(dynsyslab_double_pendulum_analytic_tests PRIVATE dynsyslab_foundation)
  add_test(NAME dynsyslab_double_pendulum_analytic_tests COMMAND dynsyslab_double_pendulum_analytic_tests)
  add_executable(dynsyslab_nonode_analytic_tests tests/nonode_analytic_tests.cpp)
  target_link_libraries(dynsyslab_nonode_analytic_tests PRIVATE dynsyslab_delay dynsyslab_pde)
  add_test(NAME dynsyslab_nonode_analytic_tests COMMAND dynsyslab_nonode_analytic_tests)
  add_executable(dynsyslab_delay_tests tests/delay_tests.cpp)
  target_link_libraries(dynsyslab_delay_tests PRIVATE dynsyslab_delay)
  add_test(NAME dynsyslab_delay_tests COMMAND dynsyslab_delay_tests)
  add_executable(dynsyslab_lattice_tests tests/lattice_tests.cpp)
  target_link_libraries(dynsyslab_lattice_tests PRIVATE dynsyslab_lattice)
  add_test(NAME dynsyslab_lattice_tests COMMAND dynsyslab_lattice_tests)
  add_executable(dynsyslab_pde_tests tests/pde_tests.cpp)
  target_link_libraries(dynsyslab_pde_tests PRIVATE dynsyslab_pde)
  add_test(NAME dynsyslab_pde_tests COMMAND dynsyslab_pde_tests)
  if(DYNSYSLAB_ENABLE_ARCHIVE_IO)
    add_executable(dynsyslab_archive_tests tests/archive_tests.cpp)
    target_link_libraries(dynsyslab_archive_tests PRIVATE dynsyslab_operations)
    add_test(NAME dynsyslab_archive_tests COMMAND dynsyslab_archive_tests)
    add_executable(dynsyslab_legacy_csv_tests tests/legacy_csv_tests.cpp)
    target_link_libraries(dynsyslab_legacy_csv_tests PRIVATE dynsyslab_foundation)
    add_test(NAME dynsyslab_legacy_csv_tests COMMAND dynsyslab_legacy_csv_tests)
  endif()
  if(DYNSYSLAB_ENABLE_LINALG)
    if(DYNSYSLAB_ENABLE_ARCHIVE_IO)
      # Deployment client deliberately excludes the umbrella core and all fitters.
      add_executable(dynsyslab_observation_runtime_tests tests/observation_runtime_tests.cpp)
      target_link_libraries(dynsyslab_observation_runtime_tests PRIVATE dynsyslab_observation dynsyslab_control)
      add_test(NAME dynsyslab_observation_runtime_tests COMMAND dynsyslab_observation_runtime_tests)
    endif()
    foreach(test phase_space phase_surface phase_forcing phase_bifurcation observable_expression echo_jet spectral sindy reconstruction koopman_control scattered_eigenfunction
                 boundary_eigenfunction eigenvalue_optimization classical_control
                 estimation decomposition history_koopman history_input_regression dmd_control edmd_control polynomial_construction eigenfunction_contract scientific_compatibility coordinate_runtime)
      add_executable(dynsyslab_${test}_tests tests/${test}_tests.cpp)
      target_link_libraries(dynsyslab_${test}_tests PRIVATE dynsyslab_core)
      add_test(NAME dynsyslab_${test}_tests COMMAND dynsyslab_${test}_tests)
    endforeach()
    add_executable(dynsyslab_linalg_c_api_tests tests/linalg_c_api_tests.cpp)
    target_link_libraries(dynsyslab_linalg_c_api_tests PRIVATE dynsyslab_linalg)
    add_test(NAME dynsyslab_linalg_c_api_tests COMMAND dynsyslab_linalg_c_api_tests)
  endif()
  if(DYNSYSLAB_ENABLE_SQLITE)
    add_executable(dynsyslab_job_store_tests tests/job_store_tests.cpp)
    target_link_libraries(dynsyslab_job_store_tests PRIVATE dynsyslab_operations)
    add_test(NAME dynsyslab_job_store_tests COMMAND dynsyslab_job_store_tests)
    if(DYNSYSLAB_ENABLE_ARCHIVE_IO)
      add_executable(dynsyslab_archive_index_tests tests/archive_index_tests.cpp)
      target_link_libraries(dynsyslab_archive_index_tests PRIVATE dynsyslab_operations)
      add_test(NAME dynsyslab_archive_index_tests COMMAND dynsyslab_archive_index_tests)
    endif()
  endif()
  add_executable(dynsyslab_c_api_tests tests/c_api_tests.cpp)
  target_link_libraries(dynsyslab_c_api_tests PRIVATE dynsyslab_c_api)
  add_test(NAME dynsyslab_c_api_tests COMMAND dynsyslab_c_api_tests)
  add_executable(dynsyslab_training_domain_tests tests/training_domain_tests.cpp)
  target_link_libraries(dynsyslab_training_domain_tests PRIVATE dynsyslab_foundation)
  add_test(NAME dynsyslab_training_domain_tests COMMAND dynsyslab_training_domain_tests)
  if(COMMAND dynsyslab_add_extension)
    dynsyslab_add_extension(dynsyslab_test_extension NO_INSTALL
      SOURCES tests/test_extension.cpp)
  else()
    add_library(dynsyslab_test_extension SHARED tests/test_extension.cpp)
    target_include_directories(dynsyslab_test_extension PRIVATE include)
    target_compile_features(dynsyslab_test_extension PRIVATE cxx_std_17)
  endif()
  add_executable(dynsyslab_extension_tests tests/extension_tests.cpp)
  target_link_libraries(dynsyslab_extension_tests PRIVATE dynsyslab_foundation)
  add_test(NAME dynsyslab_extension_tests COMMAND dynsyslab_extension_tests
    $<TARGET_FILE:dynsyslab_test_extension>)
endif()

if(DYNSYSLAB_BUILD_PYTHON)
  find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
  find_package(pybind11 CONFIG REQUIRED)
  pybind11_add_module(_dynsyslab_native MODULE src/bindings.cpp)
  target_link_libraries(_dynsyslab_native PRIVATE dynsyslab_core)
  target_compile_features(_dynsyslab_native PRIVATE cxx_std_17)
  install(TARGETS _dynsyslab_native LIBRARY DESTINATION . RUNTIME DESTINATION .)
endif()

set(_dsl_modules dynsyslab_foundation dynsyslab_delay dynsyslab_lattice dynsyslab_pde
  dynsyslab_koopman dynsyslab_c_api dynsyslab_identification)
if(TARGET dynsyslab_eigenfunctions)
  list(APPEND _dsl_modules dynsyslab_eigenfunctions)
endif()
if(TARGET dynsyslab_phase_space)
  list(APPEND _dsl_modules dynsyslab_phase_space)
endif()
if(TARGET dynsyslab_numerics)
  list(APPEND _dsl_modules dynsyslab_numerics)
endif()
if(TARGET dynsyslab_spectral)
  list(APPEND _dsl_modules dynsyslab_spectral)
endif()
if(TARGET dynsyslab_decomposition)
  list(APPEND _dsl_modules dynsyslab_decomposition)
endif()
if(TARGET dynsyslab_history_koopman)
  list(APPEND _dsl_modules dynsyslab_history_koopman dynsyslab_observation)
endif()
if(TARGET dynsyslab_reconstruction)
  list(APPEND _dsl_modules dynsyslab_reconstruction)
endif()
if(TARGET dynsyslab_pca)
  list(APPEND _dsl_modules dynsyslab_pca)
endif()
if(TARGET dynsyslab_control)
  list(APPEND _dsl_modules dynsyslab_control)
endif()
if(TARGET dynsyslab_estimation)
  list(APPEND _dsl_modules dynsyslab_estimation)
endif()
if(TARGET dynsyslab_control_maps)
  list(APPEND _dsl_modules dynsyslab_control_maps)
endif()
if(TARGET dynsyslab_control_screening)
  list(APPEND _dsl_modules dynsyslab_control_screening)
endif()
if(TARGET dynsyslab_linalg)
  list(APPEND _dsl_modules dynsyslab_linalg)
endif()
if(TARGET dynsyslab_operations)
  list(APPEND _dsl_modules dynsyslab_operations)
endif()
install(TARGETS ${_dsl_modules} ARCHIVE DESTINATION lib LIBRARY DESTINATION bin
  RUNTIME DESTINATION bin)
install(DIRECTORY include/ DESTINATION include)
install(FILES ../cmake/DynSysLabExtensions.cmake DESTINATION lib/cmake/DynSysLab)
if(TARGET dynsyslab_native_cli)
  install(TARGETS dynsyslab_native_cli RUNTIME DESTINATION bin)
endif()
