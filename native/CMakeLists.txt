cmake_minimum_required(VERSION 3.20)
project(dynsyslab_native VERSION 0.10.1 LANGUAGES CXX)

include(CTest)
option(DYNSYSLAB_BUILD_PYTHON "Build optional Python bindings" OFF)
option(DYNSYSLAB_BUILD_CLI "Build the small core diagnostic CLI" ON)
option(DYNSYSLAB_BUILD_TESTS "Build native-core tests" ON)
option(DYNSYSLAB_BUILD_SHARED_MODULES "Build versioned native DLL modules" ON)
option(DYNSYSLAB_ENABLE_ARCHIVE_IO "Build NPZ/LZMA archive compatibility" OFF)
option(DYNSYSLAB_ENABLE_LINALG "Build Eigen-backed scientific methods" OFF)
option(DYNSYSLAB_ENABLE_SUNDIALS "Build the SUNDIALS BDF backend" OFF)
option(DYNSYSLAB_ENABLE_SQLITE "Build the durable native job store" OFF)
option(DYNSYSLAB_WITH_MPI "Build the MPI halo exchange for partitioned PDE solves" OFF)
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

# Build identity recorded in archives and data packs (one translation unit, set at configure).
find_package(Git QUIET)
set(DYNSYSLAB_GIT_COMMIT "unknown")
set(DYNSYSLAB_GIT_MODIFIED 0)
if(GIT_FOUND)
  execute_process(COMMAND "${GIT_EXECUTABLE}" rev-parse HEAD WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    OUTPUT_VARIABLE _commit OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET RESULT_VARIABLE _commit_result)
  if(_commit_result EQUAL 0)
    set(DYNSYSLAB_GIT_COMMIT "${_commit}")
    execute_process(COMMAND "${GIT_EXECUTABLE}" status --porcelain --untracked-files=no WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      OUTPUT_VARIABLE _status OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    if(_status)
      set(DYNSYSLAB_GIT_MODIFIED 1)
    endif()
  endif()
endif()
set_property(SOURCE src/BuildInfo.cpp APPEND PROPERTY COMPILE_DEFINITIONS
  "DYNSYSLAB_SOURCE_COMMIT=\"${DYNSYSLAB_GIT_COMMIT}\"" "DYNSYSLAB_SOURCE_MODIFIED=${DYNSYSLAB_GIT_MODIFIED}"
  "DYNSYSLAB_COMPILER_ID=\"${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}\"" "DYNSYSLAB_BUILD_TYPE=\"$<CONFIG>\"")

add_library(dynsyslab_foundation ${DYNSYSLAB_LIBRARY_KIND}
  src/Models.cpp src/Execution.cpp src/SignalMetrics.cpp src/Npy.cpp src/Json.cpp src/Text.cpp src/FileDigest.cpp src/BuildInfo.cpp src/PortableNormal.cpp src/OutputPath.cpp
  src/LegacyCsv.cpp src/LegacyCsvImport.cpp
  src/TimeSeriesRecord.cpp src/DataObject.cpp src/Spectrum.cpp
  src/TimeSeries.cpp src/DelaySampling.cpp src/TimeSeriesOperations.cpp src/Fourier.cpp src/DiscreteFourier.cpp src/Extrema.cpp src/SpectrumComparison.cpp src/NumericalDifferentiation.cpp src/NumericalDifferentiationSmoothers.cpp src/FiniteDifferenceStencils.cpp src/TimeSeriesArchive.cpp
  src/Extension.cpp src/ExtensionModels.cpp src/ODE.cpp src/Integrator.cpp src/ControllerCallLog.cpp src/IntegratorExplicit.cpp src/IntegratorRadau.cpp src/DimensionalSystem.cpp src/PresetMetadata.cpp src/DelayRollout.cpp src/KnownInput.cpp src/InputAffineSensor.cpp src/Hierarchy.cpp
  src/AnalyticExpression.cpp src/AnalyticExpressionText.cpp src/EllipticIntegrals.cpp src/FunctionExpansion.cpp src/ExpressionTerms.cpp src/ExpansionDictionary.cpp src/ObservableExpressions.cpp src/OscillatorObservables.cpp src/BasisFunctions.cpp
  src/SpecialFunctions.cpp src/Lauricella.cpp src/QuarticOrbit.cpp src/AnalyticSolution.cpp src/AnalyticHistories.cpp src/PeriodicAnalysis.cpp src/TorusAnalysis.cpp src/CoefficientShift.cpp
  src/AnalyticQuadrature.cpp src/AnalyticDoublePendulum.cpp src/AnalyticModelLimits.cpp src/AnalyticKuramoto.cpp
  src/NeuralKnowledge.cpp src/CoupledKnowledge.cpp src/BubbleKnowledge.cpp
  src/AnalyticHenonHeiles.cpp src/AnalyticOscillators.cpp src/AnalyticRPE.cpp src/AnalyticOtherSystems.cpp
  src/SystemsHarmonic.cpp src/ResearchPapers.cpp src/ResearchPapersExtra.cpp src/ResearchPapersOscillators.cpp src/DynamicalKnowledge.cpp src/SystemEigenfunctionLibrary.cpp src/RpeDynamicalKnowledge.cpp src/DoublePendulumKnowledge.cpp
  src/Hankel.cpp src/PathSimplification.cpp src/ObservableDerivative.cpp src/Observable.cpp src/ObservableRepresentation.cpp src/FunctionDomain.cpp src/DomainGrid.cpp src/GridFieldArchive.cpp src/RecordedArchive.cpp src/TimeGridArchive.cpp src/ArchiveCompression.cpp src/System.cpp src/SystemSupport.cpp src/SystemLabels.cpp src/Targets.cpp
  src/TrainingDomain.cpp
  src/SystemsBasic.cpp src/SystemsTriangular.cpp src/SystemsMultiwell.cpp src/SystemsIntegrable.cpp src/AnalyticIntegrable.cpp src/AnalyticCelestial.cpp src/SystemsToda.cpp src/AnalyticToda.cpp src/Toda.cpp src/SystemsMechanical.cpp src/SystemsBiological.cpp
  src/SystemsBubble.cpp src/SystemsRPE.cpp src/SystemsRPEBaseline.cpp src/SystemsMarmottant.cpp src/AcousticReceiver.cpp
  src/SystemsReddySzeri.cpp)
dynsyslab_module(dynsyslab_foundation)
find_package(boost_math CONFIG REQUIRED)
find_package(boost_multiprecision CONFIG REQUIRED)
target_link_libraries(dynsyslab_foundation PRIVATE Boost::math)

# Exact geometric analysis is independent of integration, fitting, and archives.
add_library(dynsyslab_math ${DYNSYSLAB_LIBRARY_KIND}
  src/math/Arithmetic.cpp src/math/Topology.cpp src/math/HarmonicAnalysis.cpp src/math/AnalyticFourier.cpp src/math/KeplerFourier.cpp src/math/PeriodicSnapshots.cpp src/math/InvariantProjection.cpp src/math/DiscreteCycles.cpp src/math/CyclicCocycle.cpp
  src/math/DirectScattering.cpp src/math/ReflectionlessScattering.cpp src/math/ScatteringQuery.cpp
  src/math/Resonance.cpp src/math/Distributions.cpp src/math/SpectralFilters.cpp src/math/PeriodicAtlas.cpp src/math/PeriodicChart.cpp src/math/PeriodicClock.cpp src/math/TorusFourier.cpp src/math/PhaseCoverFourier.cpp src/math/TorusClock.cpp src/math/Perturbation.cpp src/math/HamiltonianPerturbation.cpp src/math/TorusQuery.cpp src/math/Query.cpp
  src/math/OrbitIntegrator.cpp src/math/OrbitFamily.cpp src/math/DampedOrbit.cpp src/math/MelnikovTools.cpp src/math/OrbitGramian.cpp src/math/OrbitQuery.cpp
  src/math/EigenvalueLattice.cpp src/math/ShellPacket.cpp src/math/ReebGraph.cpp src/math/SpectralDiagnosticsQuery.cpp
  src/math/ResolventLinearAlgebra.cpp src/math/LaplaceTransform.cpp src/math/KoopmanResolvent.cpp src/math/GeneralizedLaplace.cpp
  src/math/KoopmanFrequencyResponse.cpp src/math/ShellResolvent.cpp src/math/FrequencyGramian.cpp src/math/ResolventQuery.cpp
  src/math/ResolventLeastSquares.cpp src/math/DiscreteResolvent.cpp src/math/PronyIdentification.cpp src/math/ForcedResponseEstimation.cpp
  src/math/IsochronDecomposition.cpp src/math/HarmonicContinuation.cpp src/math/ContractionCheck.cpp src/math/ResolventDataQuery.cpp)
dynsyslab_module(dynsyslab_math)
target_link_libraries(dynsyslab_math PUBLIC dynsyslab_foundation PRIVATE Boost::multiprecision Boost::math)

add_library(dynsyslab_lie ${DYNSYSLAB_LIBRARY_KIND}
  src/lie/Polynomial.cpp src/lie/Closure.cpp src/lie/Structure.cpp src/lie/Query.cpp src/lie/Bounds.cpp
  src/lie/InverseRadius.cpp src/lie/SinusoidalMechanics.cpp)
dynsyslab_module(dynsyslab_lie)
target_link_libraries(dynsyslab_lie PUBLIC dynsyslab_foundation PRIVATE Boost::multiprecision)
if(BUILD_TESTING AND DYNSYSLAB_BUILD_TESTS)
  add_executable(dynsyslab_known_input_tests tests/known_input_tests.cpp)
  target_link_libraries(dynsyslab_known_input_tests PRIVATE dynsyslab_foundation)
  add_test(NAME dynsyslab_known_input_tests COMMAND dynsyslab_known_input_tests)
  add_executable(dynsyslab_integrable_system_tests tests/integrable_system_tests.cpp)
  add_executable(dynsyslab_toda_tests tests/toda_tests.cpp)
  target_link_libraries(dynsyslab_toda_tests PRIVATE dynsyslab_foundation)
  add_test(NAME dynsyslab_toda_tests COMMAND dynsyslab_toda_tests)
  add_executable(dynsyslab_resolvent_tools_tests tests/resolvent_tools_tests.cpp)
  target_link_libraries(dynsyslab_resolvent_tools_tests PRIVATE dynsyslab_math)
  add_test(NAME dynsyslab_resolvent_tools_tests COMMAND dynsyslab_resolvent_tools_tests)
  add_executable(dynsyslab_resolvent_data_tests tests/resolvent_data_tests.cpp)
  target_link_libraries(dynsyslab_resolvent_data_tests PRIVATE dynsyslab_math)
  add_test(NAME dynsyslab_resolvent_data_tests COMMAND dynsyslab_resolvent_data_tests)
  add_executable(dynsyslab_orbit_tools_tests tests/orbit_tools_tests.cpp)
  target_link_libraries(dynsyslab_orbit_tools_tests PRIVATE dynsyslab_math)
  add_test(NAME dynsyslab_orbit_tools_tests COMMAND dynsyslab_orbit_tools_tests)
  add_executable(dynsyslab_spectral_diagnostics_tests tests/spectral_diagnostics_tests.cpp)
  target_link_libraries(dynsyslab_spectral_diagnostics_tests PRIVATE dynsyslab_math)
  add_test(NAME dynsyslab_spectral_diagnostics_tests COMMAND dynsyslab_spectral_diagnostics_tests)
  add_executable(dynsyslab_torus_analysis_tests tests/torus_analysis_tests.cpp)
  add_executable(dynsyslab_inverse_scattering_tests tests/inverse_scattering_tests.cpp)
  target_link_libraries(dynsyslab_inverse_scattering_tests PRIVATE dynsyslab_math)
  add_test(NAME dynsyslab_inverse_scattering_tests COMMAND dynsyslab_inverse_scattering_tests)
  target_link_libraries(dynsyslab_torus_analysis_tests PRIVATE dynsyslab_math)
  add_test(NAME dynsyslab_torus_analysis_tests COMMAND dynsyslab_torus_analysis_tests)
  target_link_libraries(dynsyslab_integrable_system_tests PRIVATE dynsyslab_foundation)
  add_test(NAME dynsyslab_integrable_system_tests COMMAND dynsyslab_integrable_system_tests)
  add_executable(dynsyslab_math_tests tests/math_tests.cpp)
  add_executable(dynsyslab_torus_clock_tests tests/torus_clock_tests.cpp)
  target_link_libraries(dynsyslab_torus_clock_tests PRIVATE dynsyslab_math)
  add_test(NAME dynsyslab_torus_clock_tests COMMAND dynsyslab_torus_clock_tests)
  add_executable(dynsyslab_discrete_cycle_tests tests/discrete_cycle_tests.cpp)
  target_link_libraries(dynsyslab_discrete_cycle_tests PRIVATE dynsyslab_math)
  add_test(NAME dynsyslab_discrete_cycle_tests COMMAND dynsyslab_discrete_cycle_tests)
  target_link_libraries(dynsyslab_math_tests PRIVATE dynsyslab_math)
  add_test(NAME dynsyslab_math_tests COMMAND dynsyslab_math_tests)
  add_executable(dynsyslab_spectral_filter_tests tests/spectral_filter_tests.cpp)
  add_executable(dynsyslab_periodic_clock_tests tests/periodic_clock_tests.cpp)
  target_link_libraries(dynsyslab_periodic_clock_tests PRIVATE dynsyslab_math)
  add_test(NAME dynsyslab_periodic_clock_tests COMMAND dynsyslab_periodic_clock_tests)
  target_link_libraries(dynsyslab_spectral_filter_tests PRIVATE dynsyslab_math)
  add_test(NAME dynsyslab_spectral_filter_tests COMMAND dynsyslab_spectral_filter_tests)
  add_executable(dynsyslab_multiwell_system_tests tests/multiwell_system_tests.cpp)
  target_link_libraries(dynsyslab_multiwell_system_tests PRIVATE dynsyslab_foundation)
  add_test(NAME dynsyslab_multiwell_system_tests COMMAND dynsyslab_multiwell_system_tests)
  add_executable(dynsyslab_lie_tests tests/lie_tests.cpp)
  target_link_libraries(dynsyslab_lie_tests PRIVATE dynsyslab_lie)
  add_test(NAME dynsyslab_lie_tests COMMAND dynsyslab_lie_tests)
  add_executable(dynsyslab_sinusoidal_closure_tests tests/sinusoidal_closure_tests.cpp)
  target_link_libraries(dynsyslab_sinusoidal_closure_tests PRIVATE dynsyslab_lie)
  add_test(NAME dynsyslab_sinusoidal_closure_tests COMMAND dynsyslab_sinusoidal_closure_tests)
endif()

add_library(dynsyslab_delay ${DYNSYSLAB_LIBRARY_KIND} src/Delay.cpp
  src/AnalyticDelay.cpp src/DelayExponential.cpp src/DelayResearch.cpp)
dynsyslab_module(dynsyslab_delay)
target_link_libraries(dynsyslab_delay PUBLIC dynsyslab_foundation PRIVATE Boost::math Boost::multiprecision)

add_library(dynsyslab_lattice ${DYNSYSLAB_LIBRARY_KIND} src/Lattice.cpp src/LatticeResearch.cpp)
dynsyslab_module(dynsyslab_lattice)
target_link_libraries(dynsyslab_lattice PUBLIC dynsyslab_foundation)

add_library(dynsyslab_pde ${DYNSYSLAB_LIBRARY_KIND} src/PDE.cpp src/PDEAdaptive.cpp src/PDEParallel.cpp src/AnalyticHeat.cpp src/HeatResearch.cpp
  src/PDEHaloExchange.cpp src/PDEAdaptivePartition.cpp)
dynsyslab_module(dynsyslab_pde)
target_link_libraries(dynsyslab_pde PUBLIC dynsyslab_foundation)
if(DYNSYSLAB_WITH_MPI)
  find_package(MPI REQUIRED COMPONENTS CXX)
  target_sources(dynsyslab_pde PRIVATE src/PDEHaloExchangeMpi.cpp)
  target_link_libraries(dynsyslab_pde PUBLIC MPI::MPI_CXX)
else()
  target_sources(dynsyslab_pde PRIVATE src/PDEHaloExchangeNoMpi.cpp)
endif()

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
  src/SpectrumArchive.cpp
  src/ExpansionProjection.cpp src/FunctionInterpolation.cpp src/ScatteredInterpolants.cpp src/SpectralDictionary.cpp src/SpectralKromSerialization.cpp src/EigenfunctionFit.cpp src/PolynomialBasis.cpp src/KnownEigenfunctions.cpp src/EigenfunctionLibrary.cpp
  src/PolynomialEigenfunctionConstruction.cpp)
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
    src/EigenNumericsSvd.cpp src/EigenNumericsLeastSquaresReal.cpp src/EigenNumericsStreamingQr.cpp
    src/EigenNumericsLeastSquaresComplex.cpp src/NonlinearLeastSquares.cpp
    src/EigenNumericsNormalEquations.cpp
    src/EigenNumericsNullspace.cpp src/EigenNumericsComplexSpectrum.cpp
    src/EigenNumericsEigen.cpp src/EigenNumericsControl.cpp
    src/EigenNumericsControlFinite.cpp
    src/EigenNumericsMatrixFunctions.cpp src/EigenNumericsDiscreteGramian.cpp src/EigenNumericsHermitian.cpp
    src/EigenNumericsProductReal.cpp src/EigenNumericsProductComplex.cpp
    # Classical linear-systems kernels (core/StateSpace.hpp, BalancedTruncation.hpp,
    # SystemStructure.hpp, ControlDesign.hpp, SystemNorms.hpp): docs/architecture/LINEAR_SYSTEMS_TOOLKIT.md
    src/EigenNumericsStateSpace.cpp src/EigenNumericsStateSpaceSimulation.cpp
    src/EigenNumericsBalancedTruncation.cpp src/EigenNumericsSystemStructure.cpp
    src/EigenNumericsControlDesign.cpp src/EigenNumericsSystemNorms.cpp src/LoopMargins.cpp
    src/TransferFunction.cpp)
  dynsyslab_eigen_module(dynsyslab_numerics)
  target_link_libraries(dynsyslab_math PRIVATE dynsyslab_numerics)
  target_compile_definitions(dynsyslab_math PRIVATE DYNSYSLAB_MATH_NUMERICS=1)
  if(BUILD_TESTING AND DYNSYSLAB_BUILD_TESTS)
    target_compile_definitions(dynsyslab_resolvent_tools_tests PRIVATE DYNSYSLAB_MATH_NUMERICS=1)
  endif()
  add_library(dynsyslab_control_certificates ${DYNSYSLAB_LIBRARY_KIND}
    src/DenseSdp.cpp src/BilinearKrom.cpp src/ResidualSector.cpp src/SupportedResidualSector.cpp src/CertifiedDesign.cpp
    src/math/ControlCertificateQuery.cpp)
  dynsyslab_eigen_module(dynsyslab_control_certificates)
  target_link_libraries(dynsyslab_control_certificates PUBLIC dynsyslab_foundation)
  target_link_libraries(dynsyslab_math PRIVATE dynsyslab_control_certificates)
  target_compile_definitions(dynsyslab_math PRIVATE DYNSYSLAB_MATH_CONTROL_CERTIFICATES=1)
  if(BUILD_TESTING AND DYNSYSLAB_BUILD_TESTS)
    add_executable(dynsyslab_control_certificate_tests tests/control_certificate_tests.cpp)
    target_link_libraries(dynsyslab_control_certificate_tests PRIVATE dynsyslab_control_certificates)
    add_test(NAME dynsyslab_control_certificate_tests COMMAND dynsyslab_control_certificate_tests)
  endif()

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
    src/SINDy.cpp src/SindyLibrary.cpp src/SindyEigenvalues.cpp src/SindyEigenfunctionFit.cpp src/SindyDynamics.cpp src/SindyDynamicsIntegration.cpp src/ForecastScore.cpp src/DynamicsLibrary.cpp src/ModelLibraries.cpp src/RunInvariant.cpp
    src/BoundaryEigenfunctions.cpp src/EigenvalueOptimization.cpp
    src/FieldModelCore.cpp src/CompactRadialField.cpp src/LocalTangentField.cpp src/PairedPhaseField.cpp
    src/StructuredFields.cpp src/WaveformFit.cpp src/FieldCalibration.cpp)
  dynsyslab_module(dynsyslab_identification)
  target_link_libraries(dynsyslab_identification PUBLIC dynsyslab_eigenfunctions
    PRIVATE dynsyslab_numerics)

  # Lagrangian coherent structures and spectral submanifolds: docs/architecture/LCS_AND_SSM.md
  add_library(dynsyslab_invariant_geometry ${DYNSYSLAB_LIBRARY_KIND}
    src/FiniteTimeLyapunov.cpp src/LcsExtraction.cpp src/PowerSeries.cpp src/SsmParameterization.cpp src/SsmLearning.cpp)
  dynsyslab_module(dynsyslab_invariant_geometry)
  target_link_libraries(dynsyslab_invariant_geometry PUBLIC dynsyslab_foundation
    PRIVATE dynsyslab_eigenfunctions dynsyslab_numerics)
else()
  add_library(dynsyslab_identification ${DYNSYSLAB_LIBRARY_KIND}
    src/SINDyUnavailable.cpp src/SindyEigenfunctionFit.cpp)
  dynsyslab_module(dynsyslab_identification)
  target_link_libraries(dynsyslab_identification PUBLIC dynsyslab_eigenfunctions)
endif()

add_library(dynsyslab_koopman ${DYNSYSLAB_LIBRARY_KIND}
  src/KROM.cpp src/KoopmanModel.cpp src/KoopmanModelConstruction.cpp src/KoopmanValidation.cpp src/KoopmanDifferential.cpp
  src/KoopmanDatasets.cpp src/KoopmanInspection.cpp
  src/KoopmanModelEigenfunctions.cpp src/HankelRepresentations.cpp
  src/ScientificCompatibility.cpp)
dynsyslab_module(dynsyslab_koopman)
target_link_libraries(dynsyslab_koopman PUBLIC dynsyslab_foundation dynsyslab_eigenfunctions)
if(NOT DYNSYSLAB_ENABLE_LINALG)
  target_sources(dynsyslab_koopman PRIVATE src/SpectralUnavailable.cpp)
endif()
if(DYNSYSLAB_ENABLE_ARCHIVE_IO)
  target_sources(dynsyslab_koopman PRIVATE src/KoopmanArchive.cpp src/KromArchiveParts.cpp src/HankelBasis.cpp
    src/ControlMapArchive.cpp src/KoopmanInputMaps.cpp src/KoopmanBinding.cpp
    src/SpectralKromArchive.cpp)
endif()

if(DYNSYSLAB_ENABLE_LINALG)
  add_library(dynsyslab_spectral ${DYNSYSLAB_LIBRARY_KIND}
    src/Spectral.cpp src/DmdControl.cpp src/DmdControlEigenfunctions.cpp src/EdmdControl.cpp
    src/KoopmanSpectralFit.cpp src/KoopmanResiduals.cpp src/KoopmanMeasures.cpp src/KoopmanStochastic.cpp src/KoopmanRkhs.cpp
      src/KoopmanDictionaryFit.cpp src/PeriodicSpectralComparison.cpp src/PeriodicDelayForecast.cpp src/PeriodicPacketFeatures.cpp)
  dynsyslab_module(dynsyslab_spectral)
  target_link_libraries(dynsyslab_spectral PUBLIC
    dynsyslab_foundation dynsyslab_koopman PRIVATE dynsyslab_numerics dynsyslab_math)
  if(DYNSYSLAB_ENABLE_ARCHIVE_IO)
    target_sources(dynsyslab_spectral PRIVATE src/EdmdControlArchive.cpp)
  endif()

  add_library(dynsyslab_decomposition ${DYNSYSLAB_LIBRARY_KIND}
    src/Decomposition.cpp)
  dynsyslab_module(dynsyslab_decomposition)
  target_link_libraries(dynsyslab_decomposition PUBLIC dynsyslab_foundation
    PRIVATE dynsyslab_numerics)

  add_library(dynsyslab_observation ${DYNSYSLAB_LIBRARY_KIND}
    src/HankelOperator.cpp src/PathCrossings.cpp src/HistoryProjection.cpp src/FrozenHistoryTransform.cpp src/CoordinateTransform.cpp src/CoordinateMaps.cpp src/OrbitGeometry.cpp
    src/ReceiverArray.cpp src/NeighborReadout.cpp src/TangentReadout.cpp src/DelayDiagnostics.cpp src/AffineDynamics.cpp src/TrajectoryAudit.cpp
    src/CoordinateClosure.cpp src/jet/GlobalJetClosure.cpp src/jet/JetDecision.cpp)
  dynsyslab_eigen_module(dynsyslab_observation)
  target_link_libraries(dynsyslab_observation PUBLIC dynsyslab_foundation
    PRIVATE dynsyslab_numerics Boost::multiprecision ${CMAKE_DL_LIBS})
  if(DYNSYSLAB_ENABLE_ARCHIVE_IO)
    target_sources(dynsyslab_observation PRIVATE src/HistoryObservationArchive.cpp)
  endif()

  add_library(dynsyslab_history_koopman ${DYNSYSLAB_LIBRARY_KIND}
    src/HistoryKoopman.cpp src/HistoryInputRegression.cpp)
  dynsyslab_eigen_module(dynsyslab_history_koopman)
  target_link_libraries(dynsyslab_history_koopman PUBLIC
    dynsyslab_koopman dynsyslab_observation PRIVATE dynsyslab_numerics)

  add_library(dynsyslab_reconstruction ${DYNSYSLAB_LIBRARY_KIND}
    src/Reconstruction.cpp src/ReconstructionQuality.cpp)
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
    src/ControlRegression.cpp src/SymmetricPulseControlMap.cpp src/PathSampling.cpp src/ControlMapAudit.cpp)
  dynsyslab_module(dynsyslab_control_maps)
  target_link_libraries(dynsyslab_control_maps PUBLIC dynsyslab_koopman
    PRIVATE dynsyslab_numerics)

  add_library(dynsyslab_control_screening ${DYNSYSLAB_LIBRARY_KIND}
    src/KromScreening.cpp)
  dynsyslab_module(dynsyslab_control_screening)
  target_link_libraries(dynsyslab_control_screening PUBLIC
    dynsyslab_koopman dynsyslab_spectral)

  add_library(dynsyslab_linear_analysis ${DYNSYSLAB_LIBRARY_KIND}
    src/Observability.cpp src/ObservabilityParticipation.cpp src/ModalObservability.cpp src/Pbh.cpp src/Controllability.cpp src/KoopmanLinearAnalysis.cpp
    src/SystemStructure.cpp src/KoopmanStateSpace.cpp src/EmpiricalObservability.cpp
    src/EigenfunctionObservability.cpp src/KromDifferentialAnalysis.cpp src/KromLocalControlAnalysis.cpp)
  dynsyslab_module(dynsyslab_linear_analysis)
  # The toolkit's dense kernels live in dynsyslab_numerics; their core headers are
  # part of this module's API, so numerics is a public dependency.
  target_link_libraries(dynsyslab_linear_analysis PUBLIC dynsyslab_koopman dynsyslab_numerics)

  # Feedback design and runtime: estimation (KLQE) and control (KLQR, classical)
  # in one module. Gains come from ControlDesign in dynsyslab_linear_analysis;
  # the estimator/controller split lives in the workflow job format and factories.
  add_library(dynsyslab_control ${DYNSYSLAB_LIBRARY_KIND}
    src/Estimation.cpp src/EstimationRuntime.cpp src/EstimationSelection.cpp src/EstimationRealization.cpp
    src/ReducedOrderEstimation.cpp src/EnsembleBayes.cpp src/PeriodicEstimation.cpp src/PeriodicIdentification.cpp src/PeriodicFeatureIdentification.cpp src/PeriodicTransport.cpp src/ShootingIdentification.cpp src/SequentialShooting.cpp src/OnlineShooting.cpp
    src/LQR.cpp src/IoDelayModel.cpp src/IoDelayTracker.cpp src/ClassicalControl.cpp src/ReducedOrderControl.cpp src/KoopmanBilinearBridge.cpp
    src/LinearObserver.cpp src/PodGalerkin.cpp src/EigenfunctionRefit.cpp src/ControlMapPath.cpp)
  dynsyslab_module(dynsyslab_control)
  target_link_libraries(dynsyslab_control PUBLIC
    dynsyslab_koopman dynsyslab_linear_analysis dynsyslab_math dynsyslab_control_certificates PRIVATE dynsyslab_numerics)
  if(BUILD_TESTING AND DYNSYSLAB_BUILD_TESTS)
    add_executable(dynsyslab_control_estimation_core_tests tests/control_estimation_core_tests.cpp)
    target_link_libraries(dynsyslab_control_estimation_core_tests PRIVATE dynsyslab_control
      dynsyslab_control_maps dynsyslab_reconstruction)
    add_test(NAME dynsyslab_control_estimation_core_tests COMMAND dynsyslab_control_estimation_core_tests)
    add_executable(dynsyslab_io_delay_model_tests tests/io_delay_model_tests.cpp)
    target_link_libraries(dynsyslab_io_delay_model_tests PRIVATE dynsyslab_control)
    add_test(NAME dynsyslab_io_delay_model_tests COMMAND dynsyslab_io_delay_model_tests)
    add_executable(dynsyslab_periodic_transport_tests tests/periodic_transport_tests.cpp)
    target_link_libraries(dynsyslab_periodic_transport_tests PRIVATE dynsyslab_control)
    add_test(NAME dynsyslab_periodic_transport_tests COMMAND dynsyslab_periodic_transport_tests)
    add_executable(dynsyslab_shooting_identification_tests tests/shooting_identification_tests.cpp)
    target_link_libraries(dynsyslab_shooting_identification_tests PRIVATE dynsyslab_control)
    add_test(NAME dynsyslab_shooting_identification_tests COMMAND dynsyslab_shooting_identification_tests)
    add_executable(dynsyslab_sequential_shooting_tests tests/sequential_shooting_tests.cpp)
    target_link_libraries(dynsyslab_sequential_shooting_tests PRIVATE dynsyslab_control)
    add_test(NAME dynsyslab_sequential_shooting_tests COMMAND dynsyslab_sequential_shooting_tests)
    add_executable(dynsyslab_online_shooting_tests tests/online_shooting_tests.cpp)
    target_link_libraries(dynsyslab_online_shooting_tests PRIVATE dynsyslab_control)
    add_test(NAME dynsyslab_online_shooting_tests COMMAND dynsyslab_online_shooting_tests)
    add_executable(dynsyslab_periodic_estimation_tests tests/periodic_estimation_tests.cpp)
    add_executable(dynsyslab_periodic_feature_tests tests/periodic_feature_tests.cpp)
    target_link_libraries(dynsyslab_periodic_feature_tests PRIVATE dynsyslab_control dynsyslab_spectral dynsyslab_math dynsyslab_numerics)
    add_test(NAME dynsyslab_periodic_feature_tests COMMAND dynsyslab_periodic_feature_tests)
    target_link_libraries(dynsyslab_periodic_estimation_tests PRIVATE dynsyslab_control)
    add_test(NAME dynsyslab_periodic_estimation_tests COMMAND dynsyslab_periodic_estimation_tests)
  endif()
endif()

if(DYNSYSLAB_ENABLE_SQLITE)
  find_package(unofficial-sqlite3 CONFIG REQUIRED)
  add_library(dynsyslab_operations ${DYNSYSLAB_LIBRARY_KIND}
    src/JobStore.cpp src/JobScheduling.cpp src/CampaignRegistry.cpp src/ArchiveAudit.cpp src/ArchiveDeduplication.cpp
    src/ArchiveAuditReport.cpp src/ArchiveMaintenance.cpp src/StorageLocations.cpp src/CampaignCatalog.cpp
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
if(TARGET dynsyslab_invariant_geometry)
  target_link_libraries(dynsyslab_core INTERFACE dynsyslab_invariant_geometry)
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
  if(DYNSYSLAB_ENABLE_LINALG AND DYNSYSLAB_ENABLE_ARCHIVE_IO)
    add_executable(dynsyslab_streaming_spectral_tests tests/streaming_spectral_tests.cpp)
    target_link_libraries(dynsyslab_streaming_spectral_tests PRIVATE dynsyslab_spectral dynsyslab_math dynsyslab_numerics)
    add_test(NAME dynsyslab_streaming_spectral_tests COMMAND dynsyslab_streaming_spectral_tests)
    add_executable(dynsyslab_koopman_spectral_tests tests/koopman_spectral_tests.cpp)
    target_link_libraries(dynsyslab_koopman_spectral_tests PRIVATE dynsyslab_spectral dynsyslab_math dynsyslab_numerics)
    add_test(NAME dynsyslab_koopman_spectral_tests COMMAND dynsyslab_koopman_spectral_tests)
    add_executable(dynsyslab_observable_dictionary_tests tests/observable_dictionary_tests.cpp)
    target_link_libraries(dynsyslab_observable_dictionary_tests PRIVATE dynsyslab_spectral)
    add_test(NAME dynsyslab_observable_dictionary_tests COMMAND dynsyslab_observable_dictionary_tests)
  endif()
  foreach(test core system integrator analytic_solution triangular_system)
    add_executable(dynsyslab_${test}_tests tests/${test}_tests.cpp)
    target_link_libraries(dynsyslab_${test}_tests PRIVATE dynsyslab_core)
    add_test(NAME dynsyslab_${test}_tests COMMAND dynsyslab_${test}_tests)
  endforeach()
  add_executable(dynsyslab_numerical_differentiation_tests tests/numerical_differentiation_tests.cpp)
  target_link_libraries(dynsyslab_numerical_differentiation_tests PRIVATE dynsyslab_foundation)
  add_test(NAME dynsyslab_numerical_differentiation_tests COMMAND dynsyslab_numerical_differentiation_tests)
  set_tests_properties(dynsyslab_numerical_differentiation_tests PROPERTIES LABELS "native;differentiation;noise")
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
  add_executable(dynsyslab_integrator_methods_tests tests/integrator_methods_tests.cpp)
  target_link_libraries(dynsyslab_integrator_methods_tests PRIVATE dynsyslab_pde dynsyslab_delay dynsyslab_lattice)
  add_test(NAME dynsyslab_integrator_methods_tests COMMAND dynsyslab_integrator_methods_tests)
  add_executable(dynsyslab_pde_tests tests/pde_tests.cpp)
  target_link_libraries(dynsyslab_pde_tests PRIVATE dynsyslab_pde)
  add_test(NAME dynsyslab_pde_tests COMMAND dynsyslab_pde_tests)
  add_executable(dynsyslab_pde_parallel_tests tests/pde_parallel_tests.cpp)
  target_link_libraries(dynsyslab_pde_parallel_tests PRIVATE dynsyslab_pde)
  add_test(NAME dynsyslab_pde_parallel_tests COMMAND dynsyslab_pde_parallel_tests)
  if(DYNSYSLAB_ENABLE_ARCHIVE_IO)
    add_executable(dynsyslab_output_path_tests tests/output_path_tests.cpp)
    target_link_libraries(dynsyslab_output_path_tests PRIVATE dynsyslab_operations)
    add_test(NAME dynsyslab_output_path_tests COMMAND dynsyslab_output_path_tests)
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
    foreach(test phase_space phase_surface phase_forcing phase_bifurcation observable_expression echo_jet coordinate_closure global_jet spectral sindy reconstruction koopman_control scattered_eigenfunction eigenfunction_objects eigenfunction_expansion eigenfunction_observable system_eigenfunction_library coordinate_maps time_series_operations
                 boundary_eigenfunction eigenvalue_optimization classical_control
                 estimation estimation_runtime linear_analysis krom_audit linear_systems loop_robustness empirical_observability reduced_order decomposition history_koopman history_input_regression dmd_control edmd_control polynomial_construction eigenfunction_contract scientific_compatibility coordinate_runtime identification_core
                 lagrangian_coherence spectral_submanifold extrema orbit_geometry trajectory_kernels model_libraries field_models)
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
      add_executable(dynsyslab_archive_maintenance_tests tests/archive_maintenance_tests.cpp)
      target_link_libraries(dynsyslab_archive_maintenance_tests PRIVATE dynsyslab_operations)
      add_test(NAME dynsyslab_archive_maintenance_tests COMMAND dynsyslab_archive_maintenance_tests)
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
  # The DLL is loaded at runtime, so linking does not establish this dependency.
  add_dependencies(dynsyslab_extension_tests dynsyslab_test_extension)
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
  dynsyslab_koopman dynsyslab_c_api dynsyslab_identification dynsyslab_lie dynsyslab_math)
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
if(TARGET dynsyslab_linear_analysis)
  list(APPEND _dsl_modules dynsyslab_linear_analysis)
endif()
install(TARGETS ${_dsl_modules} ARCHIVE DESTINATION lib LIBRARY DESTINATION bin
  RUNTIME DESTINATION bin)
install(DIRECTORY include/ DESTINATION include)
install(FILES ../cmake/DynSysLabExtensions.cmake DESTINATION lib/cmake/DynSysLab)
if(TARGET dynsyslab_native_cli)
  install(TARGETS dynsyslab_native_cli RUNTIME DESTINATION bin)
endif()
