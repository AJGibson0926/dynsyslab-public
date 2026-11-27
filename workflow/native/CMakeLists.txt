cmake_minimum_required(VERSION 3.20)
project(dynsyslab_workflow VERSION 0.10.0 LANGUAGES CXX)

if(MSVC)
  set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
endif()
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

add_library(dynsyslab_workflow SHARED
  src/Run.cpp src/Plan.cpp src/SystemResolution.cpp src/SolveDispatch.cpp src/SolveSet.cpp src/Solve.cpp src/SolveSettings.cpp src/SolveCheckpoint.cpp src/SolveArchive.cpp src/SolveContinuation.cpp src/ExtensionModelSolve.cpp src/InitialConditions.cpp src/ControlFactory.cpp src/ControlFactoryKoopman.cpp src/ControlTargets.cpp src/PrescribedControl.cpp src/PrescribedInputs.cpp src/DelaySolve.cpp src/LatticeSolve.cpp src/PDESolve.cpp src/Hierarchy.cpp src/Campaign.cpp src/Experiment.cpp src/ExperimentExpansion.cpp src/Postprocess.cpp src/SindyEigenfunctions.cpp src/PipelineFingerprint.cpp
  src/WorkflowSupport.cpp src/WorkflowExtensions.cpp src/PostprocessSupport.cpp
  src/SourceCoordinates.cpp src/ObservableQuery.cpp src/CoordinateClosureQuery.cpp
  src/AnalyticQuery.cpp src/PhaseCoverQuery.cpp src/AnalyticHistoryQuery.cpp src/TorusQuery.cpp src/PeriodicQuery.cpp src/PeriodicControlQuery.cpp src/PeriodicSnapshotQuery.cpp src/PeriodicSpectralQuery.cpp src/AnalyticArchive.cpp src/LieAlgebraQuery.cpp src/KromAuditQuery.cpp src/KromTrajectoryAudit.cpp
  src/PhaseSpaceQuery.cpp src/PhaseSpaceSerialization.cpp src/PhaseForcingConfiguration.cpp
  src/KnowledgeSerialization.cpp src/ResearchCatalogue.cpp
  src/LinearAnalysisPostprocess.cpp src/LinearSystemsPostprocess.cpp src/LinearSystemsClassical.cpp src/EmpiricalObservabilityPostprocess.cpp src/ReconstructionPostprocess.cpp src/HankelAnalysis.cpp src/HankelArchive.cpp src/HdmdPostprocess.cpp src/HdmdObjects.cpp src/HdmdPooled.cpp src/CoordinateDelayTaps.cpp src/SampledDomain.cpp src/SindyTraining.cpp src/ObjectPublishing.cpp src/TimeSeriesSource.cpp src/HankelSvdPostprocess.cpp src/MeasurementMap.cpp src/ExpansionArchive.cpp src/EigenfunctionExpansionPostprocess.cpp src/DictionaryTerms.cpp src/DomainFieldPostprocess.cpp src/LcsPostprocess.cpp src/SsmPostprocess.cpp src/OrbitGeometryPostprocess.cpp src/ExperimentSweep.cpp src/ReceiverPostprocess.cpp src/TrajectoryKernelsPostprocess.cpp src/DelayDiagnosticsPostprocess.cpp src/RunInvariantPostprocess.cpp src/FieldModelPostprocess.cpp src/PulseControlPostprocess.cpp src/AnimationManifest.cpp src/EigenfunctionErrorsPostprocess.cpp src/DmdPostprocess.cpp src/DmdControlPostprocess.cpp
  src/KromControlSearch.cpp src/KromControlSearchSupport.cpp src/KromControlSearchVariants.cpp
  src/TrainingPreflight.cpp
  src/EchoHistoryKoopman.cpp src/ObservationSource.cpp src/HistoryControl.cpp src/MeasuredJetControl.cpp src/IoDelayControl.cpp src/IoDelayPostprocess.cpp src/HistoryControlRegression.cpp src/InputMap.cpp src/ModeSelection.cpp src/ModeSelectionData.cpp src/ModeSelectionScore.cpp src/EigenfunctionPolynomialFit.cpp src/KlqrAugmentation.cpp src/ControlTargetLift.cpp src/DerivativeBenchmark.cpp src/HistoryPostprocess.cpp src/CoordinateSindy.cpp
  src/EstimationPostprocess.cpp src/DecompositionPostprocess.cpp
  src/CoordinateTransformPostprocess.cpp src/CoordinateChannels.cpp src/CoordinateChangeKinds.cpp src/CoordinateResample.cpp src/CoordinateDiagnostics.cpp src/CoordinateInverse.cpp src/HankelShorthand.cpp src/SindyDynamicsPostprocess.cpp
  src/PcaPostprocess.cpp src/KromPostprocess.cpp src/KromEvidence.cpp
  src/KoopmanSpectralQuery.cpp src/KoopmanSpectralPostprocess.cpp)
target_sources(dynsyslab_workflow PRIVATE src/KromValidationPostprocess.cpp src/PeriodicIdentificationQuery.cpp src/PeriodicFeatureQuery.cpp src/PeriodicFeedback.cpp
  src/KoopmanFeedback.cpp src/EstimatedControl.cpp src/EstimatorFactory.cpp src/KlqeProcessCovariance.cpp
  src/HeldInputRecords.cpp src/BilinearCertificatePostprocess.cpp src/BilinearCertificateBinding.cpp
  src/KromDatasetSupport.cpp src/KromValidationPlan.cpp
  src/JobGroup.cpp src/JobGroupGraph.cpp src/ArchiveDigest.cpp src/PngCompare.cpp src/ViewerFigureJob.cpp src/ViewerFigureSources.cpp
  src/DataPack.cpp src/DataPackLineage.cpp
  src/KromEditing.cpp src/ReferenceSignal.cpp src/ControlScore.cpp src/ControlScoreEstimation.cpp src/FigureDataset.cpp src/FigureDatasetOperations.cpp src/BifurcationPoints.cpp src/ObservableSpectrumComparison.cpp
  src/BoundaryPostprocess.cpp src/EigenvalueOptimizationPostprocess.cpp
  src/CampaignDesign.cpp src/HistoryPolicy.cpp src/JobTemplates.cpp src/PostprocessGeometry.cpp
  src/ClassicalEstimator.cpp src/KromDomainPolicy.cpp src/ControlFactoryClassical.cpp
  src/KromModeRules.cpp src/KromValidationWindows.cpp
  src/SolveHistory.cpp
  src/InputMapPulses.cpp src/KromTools.cpp src/KlqeDesignReuse.cpp
  src/LinearAnalysisModal.cpp src/KromTrainingPlan.cpp)
target_include_directories(dynsyslab_workflow PUBLIC include)
target_link_libraries(dynsyslab_workflow PUBLIC dynsyslab_core dynsyslab_delay
  dynsyslab_lattice dynsyslab_pde PRIVATE dynsyslab_operations dynsyslab_lie dynsyslab_math)
target_compile_features(dynsyslab_workflow PUBLIC cxx_std_17)
if(WIN32)
  target_link_libraries(dynsyslab_workflow PRIVATE windowscodecs ole32)
  set_target_properties(dynsyslab_workflow PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
  target_compile_definitions(dynsyslab_workflow PUBLIC WIN32_LEAN_AND_MEAN NOMINMAX)
endif()
if(MSVC)
  target_compile_options(dynsyslab_workflow PRIVATE /W4 /permissive- /EHsc /utf-8 /MP4 /bigobj)
else()
  target_compile_options(dynsyslab_workflow PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_library(dynsyslab_queue SHARED src/Queue.cpp src/QueueSubmit.cpp src/QueueServe.cpp src/QueueContinuation.cpp src/QueueStaging.cpp
  src/CampaignStatus.cpp src/CampaignRestructure.cpp src/Benchmark.cpp)
target_include_directories(dynsyslab_queue PUBLIC include)
target_link_libraries(dynsyslab_queue PUBLIC dynsyslab_workflow dynsyslab_operations)
target_compile_features(dynsyslab_queue PUBLIC cxx_std_17)
if(WIN32)
  set_target_properties(dynsyslab_queue PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
  target_compile_definitions(dynsyslab_queue PUBLIC WIN32_LEAN_AND_MEAN NOMINMAX)
endif()
if(MSVC)
  target_compile_options(dynsyslab_queue PRIVATE /W4 /permissive- /EHsc /utf-8 /MP4)
else()
  target_compile_options(dynsyslab_queue PRIVATE -Wall -Wextra -Wpedantic)
endif()

if(BUILD_TESTING AND DYNSYSLAB_BUILD_TESTS)
  add_executable(dynsyslab_lie_query_tests tests/lie_query_tests.cpp)
  target_link_libraries(dynsyslab_lie_query_tests PRIVATE dynsyslab_workflow dynsyslab_workflow_c_api)
  get_filename_component(_dsl_lie_examples "${CMAKE_CURRENT_SOURCE_DIR}/../../examples/native/lie-algebra" ABSOLUTE)
  target_compile_definitions(dynsyslab_lie_query_tests PRIVATE DYNSYSLAB_LIE_EXAMPLES="${_dsl_lie_examples}")
  add_test(NAME dynsyslab_lie_query_tests COMMAND dynsyslab_lie_query_tests)
  add_executable(dynsyslab_sindy_dynamics_tests tests/sindy_dynamics_tests.cpp)
  target_link_libraries(dynsyslab_sindy_dynamics_tests PRIVATE dynsyslab_workflow dynsyslab_operations)
  add_test(NAME dynsyslab_sindy_dynamics_tests COMMAND dynsyslab_sindy_dynamics_tests)
  add_executable(dynsyslab_koopman_feedback_tests tests/koopman_feedback_tests.cpp)
  target_link_libraries(dynsyslab_koopman_feedback_tests PRIVATE dynsyslab_workflow)
  add_test(NAME dynsyslab_koopman_feedback_tests COMMAND dynsyslab_koopman_feedback_tests)
  add_executable(dynsyslab_pipeline_boundary_tests tests/pipeline_boundary_tests.cpp)
  target_link_libraries(dynsyslab_pipeline_boundary_tests PRIVATE dynsyslab_workflow dynsyslab_workflow_c_api)
  add_test(NAME dynsyslab_pipeline_boundary_tests COMMAND dynsyslab_pipeline_boundary_tests)
  add_executable(dynsyslab_experiment_job_tests tests/experiment_job_tests.cpp)
  target_link_libraries(dynsyslab_experiment_job_tests PRIVATE dynsyslab_workflow)
  add_test(NAME dynsyslab_experiment_job_tests COMMAND dynsyslab_experiment_job_tests)
  add_executable(dynsyslab_dmd_control_postprocess_tests tests/dmd_control_postprocess_tests.cpp)
  target_link_libraries(dynsyslab_dmd_control_postprocess_tests PRIVATE dynsyslab_workflow dynsyslab_workflow_c_api)
  add_test(NAME dynsyslab_dmd_control_postprocess_tests COMMAND dynsyslab_dmd_control_postprocess_tests)
  add_executable(dynsyslab_research_catalogue_tests tests/research_catalogue_tests.cpp)
  target_link_libraries(dynsyslab_research_catalogue_tests PRIVATE dynsyslab_workflow)
  add_test(NAME dynsyslab_research_catalogue_tests COMMAND dynsyslab_research_catalogue_tests)
  add_executable(dynsyslab_figure_pipeline_tests tests/figure_pipeline_tests.cpp)
  target_link_libraries(dynsyslab_figure_pipeline_tests PRIVATE dynsyslab_workflow)
  add_test(NAME dynsyslab_figure_pipeline_tests COMMAND dynsyslab_figure_pipeline_tests)
  add_executable(dynsyslab_job_group_tests tests/job_group_tests.cpp)
  target_link_libraries(dynsyslab_job_group_tests PRIVATE dynsyslab_workflow)
  add_test(NAME dynsyslab_job_group_tests COMMAND dynsyslab_job_group_tests)
  add_executable(dynsyslab_hankel_observables_tests tests/hankel_observables_tests.cpp)
  target_link_libraries(dynsyslab_hankel_observables_tests PRIVATE dynsyslab_workflow)
  add_test(NAME dynsyslab_hankel_observables_tests COMMAND dynsyslab_hankel_observables_tests)
  add_executable(dynsyslab_observable_pipeline_tests tests/observable_pipeline_tests.cpp)
  target_link_libraries(dynsyslab_observable_pipeline_tests PRIVATE dynsyslab_workflow)
  add_test(NAME dynsyslab_observable_pipeline_tests COMMAND dynsyslab_observable_pipeline_tests)
  add_executable(dynsyslab_sindy_inputs_tests tests/sindy_inputs_tests.cpp)
  target_link_libraries(dynsyslab_sindy_inputs_tests PRIVATE dynsyslab_workflow)
  add_test(NAME dynsyslab_sindy_inputs_tests COMMAND dynsyslab_sindy_inputs_tests)
  add_executable(dynsyslab_resdmd_krom_tests tests/resdmd_krom_tests.cpp)
  target_link_libraries(dynsyslab_resdmd_krom_tests PRIVATE dynsyslab_workflow)
  add_test(NAME dynsyslab_resdmd_krom_tests COMMAND dynsyslab_resdmd_krom_tests)
  add_executable(dynsyslab_identification_workflow_tests tests/identification_workflow_tests.cpp)
  target_link_libraries(dynsyslab_identification_workflow_tests PRIVATE dynsyslab_workflow)
  add_test(NAME dynsyslab_identification_workflow_tests COMMAND dynsyslab_identification_workflow_tests)
  add_executable(dynsyslab_coordinate_closure_query_tests
    tests/coordinate_closure_query_tests.cpp)
  target_link_libraries(dynsyslab_coordinate_closure_query_tests PRIVATE dynsyslab_workflow)
  add_test(NAME dynsyslab_coordinate_closure_query_tests
    COMMAND dynsyslab_coordinate_closure_query_tests)
  set_tests_properties(dynsyslab_coordinate_closure_query_tests PROPERTIES
    LABELS "native;workflow;observability;closure")
  add_executable(dynsyslab_solve_features_tests tests/solve_features_tests.cpp)
  target_link_libraries(dynsyslab_solve_features_tests PRIVATE dynsyslab_workflow)
  add_test(NAME dynsyslab_solve_features_tests COMMAND dynsyslab_solve_features_tests)
  set_tests_properties(dynsyslab_solve_features_tests PROPERTIES LABELS "native;workflow;solver")
  add_executable(dynsyslab_solve_history_tests tests/solve_history_tests.cpp)
  target_link_libraries(dynsyslab_solve_history_tests PRIVATE dynsyslab_workflow)
  add_test(NAME dynsyslab_solve_history_tests COMMAND dynsyslab_solve_history_tests)
  set_tests_properties(dynsyslab_solve_history_tests PROPERTIES LABELS "native;workflow;solver")
  add_executable(dynsyslab_initial_condition_tests
    tests/initial_condition_tests.cpp)
  target_link_libraries(dynsyslab_initial_condition_tests PRIVATE dynsyslab_workflow)
  add_test(NAME dynsyslab_initial_condition_tests
    COMMAND dynsyslab_initial_condition_tests)
  set_tests_properties(dynsyslab_initial_condition_tests PROPERTIES
    LABELS "native;workflow;solver;initial-conditions")
  add_executable(dynsyslab_workflow_queue_tests tests/queue_tests.cpp)
  target_link_libraries(dynsyslab_workflow_queue_tests PRIVATE dynsyslab_queue)
  add_test(NAME dynsyslab_workflow_queue_tests COMMAND dynsyslab_workflow_queue_tests)
  set_tests_properties(dynsyslab_workflow_queue_tests PROPERTIES
    LABELS "native;workflow;queue;ui" TIMEOUT 20)
  add_executable(dynsyslab_control_estimation_workflow_tests tests/control_estimation_workflow_tests.cpp)
  target_link_libraries(dynsyslab_control_estimation_workflow_tests PRIVATE dynsyslab_workflow dynsyslab_koopman)
  add_test(NAME dynsyslab_control_estimation_workflow_tests COMMAND dynsyslab_control_estimation_workflow_tests)
  set_tests_properties(dynsyslab_control_estimation_workflow_tests PROPERTIES
    LABELS "native;workflow;control;estimation;krom" TIMEOUT 300)
  add_executable(dynsyslab_rpe_klqr_pipeline_tests tests/rpe_klqr_pipeline_tests.cpp)
  target_link_libraries(dynsyslab_rpe_klqr_pipeline_tests PRIVATE dynsyslab_workflow dynsyslab_koopman)
  get_filename_component(_dsl_rpe_experiment
    "${CMAKE_CURRENT_SOURCE_DIR}/../../examples/native/experiments/rpe-viscous-klqr-multisine/base-r.json" ABSOLUTE)
  target_compile_definitions(dynsyslab_rpe_klqr_pipeline_tests PRIVATE DYNSYSLAB_RPE_EXPERIMENT="${_dsl_rpe_experiment}")
  add_test(NAME dynsyslab_rpe_klqr_pipeline_tests COMMAND dynsyslab_rpe_klqr_pipeline_tests)
  set_tests_properties(dynsyslab_rpe_klqr_pipeline_tests PROPERTIES
    LABELS "native;workflow;rpe;krom;klqr;pipeline" TIMEOUT 600)
  add_executable(dynsyslab_queue_recovery_tests tests/queue_recovery_tests.cpp)
  target_link_libraries(dynsyslab_queue_recovery_tests PRIVATE dynsyslab_queue)
  add_test(NAME dynsyslab_queue_recovery_tests COMMAND dynsyslab_queue_recovery_tests)
  set_tests_properties(dynsyslab_queue_recovery_tests PROPERTIES
    LABELS "native;workflow;queue;campaign" TIMEOUT 120)
  add_executable(dynsyslab_queue_scheduling_tests tests/queue_scheduling_tests.cpp)
  target_link_libraries(dynsyslab_queue_scheduling_tests PRIVATE dynsyslab_queue)
  add_test(NAME dynsyslab_queue_scheduling_tests COMMAND dynsyslab_queue_scheduling_tests)
  set_tests_properties(dynsyslab_queue_scheduling_tests PROPERTIES
    LABELS "native;workflow;queue;campaign" TIMEOUT 120)
  add_executable(dynsyslab_workflow_plan_tests tests/workflow_plan_tests.cpp)
  target_link_libraries(dynsyslab_workflow_plan_tests PRIVATE dynsyslab_workflow)
  add_test(NAME dynsyslab_workflow_plan_tests COMMAND dynsyslab_workflow_plan_tests)
  add_executable(dynsyslab_solve_set_tests tests/solve_set_tests.cpp)
  target_link_libraries(dynsyslab_solve_set_tests PRIVATE dynsyslab_workflow)
  add_test(NAME dynsyslab_solve_set_tests COMMAND dynsyslab_solve_set_tests)
  set_tests_properties(dynsyslab_solve_set_tests PROPERTIES
    LABELS "native;workflow;solver;ui;hierarchy")
  add_executable(dynsyslab_training_preflight_tests
    tests/training_preflight_tests.cpp)
  target_link_libraries(dynsyslab_training_preflight_tests PRIVATE dynsyslab_workflow)
  add_test(NAME dynsyslab_training_preflight_tests
    COMMAND dynsyslab_training_preflight_tests)
  set_tests_properties(dynsyslab_training_preflight_tests PROPERTIES
    LABELS "native;workflow;solver;ui;training" TIMEOUT 60)
  add_executable(dynsyslab_pca_postprocess_tests tests/pca_postprocess_tests.cpp)
  target_link_libraries(dynsyslab_pca_postprocess_tests PRIVATE dynsyslab_workflow)
  add_test(NAME dynsyslab_pca_postprocess_tests COMMAND dynsyslab_pca_postprocess_tests)
  set_tests_properties(dynsyslab_pca_postprocess_tests PROPERTIES
    LABELS "native;parity;postprocess;pca")
  add_executable(dynsyslab_krom_postprocess_tests tests/krom_postprocess_tests.cpp)
  target_link_libraries(dynsyslab_krom_postprocess_tests PRIVATE dynsyslab_workflow)
  add_test(NAME dynsyslab_krom_postprocess_tests COMMAND dynsyslab_krom_postprocess_tests)
  set_tests_properties(dynsyslab_krom_postprocess_tests PROPERTIES
    LABELS "native;workflow;postprocess;krom")
  add_executable(dynsyslab_krom_validation_tests tests/krom_validation_tests.cpp)
  target_link_libraries(dynsyslab_krom_validation_tests PRIVATE dynsyslab_workflow)
  add_test(NAME dynsyslab_krom_validation_tests COMMAND dynsyslab_krom_validation_tests)
  set_tests_properties(dynsyslab_krom_validation_tests PROPERTIES
    LABELS "native;workflow;postprocess;krom;validation")
  add_executable(dynsyslab_echo_history_koopman_tests
    tests/echo_history_koopman_tests.cpp)
  target_link_libraries(dynsyslab_echo_history_koopman_tests PRIVATE
    dynsyslab_workflow)
  add_test(NAME dynsyslab_echo_history_koopman_tests
    COMMAND dynsyslab_echo_history_koopman_tests)
  set_tests_properties(dynsyslab_echo_history_koopman_tests PROPERTIES
    LABELS "native;workflow;postprocess;koopman;echo-history")
  add_executable(dynsyslab_expansion_producers_tests tests/expansion_producers_tests.cpp)
  target_link_libraries(dynsyslab_expansion_producers_tests PRIVATE dynsyslab_workflow)
  add_test(NAME dynsyslab_expansion_producers_tests COMMAND dynsyslab_expansion_producers_tests)
  add_executable(dynsyslab_receiver_pipeline_tests tests/receiver_pipeline_tests.cpp)
  target_link_libraries(dynsyslab_receiver_pipeline_tests PRIVATE dynsyslab_workflow)
  add_test(NAME dynsyslab_receiver_pipeline_tests COMMAND dynsyslab_receiver_pipeline_tests)
  add_executable(dynsyslab_experiment_sweep_tests tests/experiment_sweep_tests.cpp)
  target_link_libraries(dynsyslab_experiment_sweep_tests PRIVATE dynsyslab_workflow)
  add_test(NAME dynsyslab_experiment_sweep_tests COMMAND dynsyslab_experiment_sweep_tests)
  add_executable(dynsyslab_geometry_postprocess_tests tests/geometry_postprocess_tests.cpp)
  target_link_libraries(dynsyslab_geometry_postprocess_tests PRIVATE dynsyslab_workflow)
  add_test(NAME dynsyslab_geometry_postprocess_tests COMMAND dynsyslab_geometry_postprocess_tests)
  add_executable(dynsyslab_field_model_workflow_tests tests/field_model_workflow_tests.cpp)
  target_link_libraries(dynsyslab_field_model_workflow_tests PRIVATE dynsyslab_workflow)
  add_test(NAME dynsyslab_field_model_workflow_tests COMMAND dynsyslab_field_model_workflow_tests)
  add_executable(dynsyslab_pulse_control_workflow_tests tests/pulse_control_workflow_tests.cpp)
  target_link_libraries(dynsyslab_pulse_control_workflow_tests PRIVATE dynsyslab_workflow)
  add_test(NAME dynsyslab_pulse_control_workflow_tests COMMAND dynsyslab_pulse_control_workflow_tests)
  add_executable(dynsyslab_domain_field_tests tests/domain_field_tests.cpp)
  target_link_libraries(dynsyslab_domain_field_tests PRIVATE dynsyslab_workflow)
  add_test(NAME dynsyslab_domain_field_tests COMMAND dynsyslab_domain_field_tests)
  add_executable(dynsyslab_coordinate_record_tests tests/coordinate_record_tests.cpp)
  target_link_libraries(dynsyslab_coordinate_record_tests PRIVATE dynsyslab_workflow)
  add_test(NAME dynsyslab_coordinate_record_tests COMMAND dynsyslab_coordinate_record_tests)
  add_executable(dynsyslab_hankel_svd_record_tests tests/hankel_svd_record_tests.cpp)
  target_link_libraries(dynsyslab_hankel_svd_record_tests PRIVATE dynsyslab_workflow)
  add_test(NAME dynsyslab_hankel_svd_record_tests COMMAND dynsyslab_hankel_svd_record_tests)
  add_executable(dynsyslab_hdmd_scattered_krom_tests
    tests/hdmd_scattered_krom_tests.cpp)
  target_link_libraries(dynsyslab_hdmd_scattered_krom_tests PRIVATE dynsyslab_workflow)
  add_test(NAME dynsyslab_hdmd_scattered_krom_tests
    COMMAND dynsyslab_hdmd_scattered_krom_tests)
  set_tests_properties(dynsyslab_hdmd_scattered_krom_tests PROPERTIES
    LABELS "native;workflow;postprocess;krom;hdmd")
  # Frozen-reference and analytic replacements for the retired MATLAB parity scripts.
  add_executable(dynsyslab_colbrook_reference_tests tests/colbrook_reference_tests.cpp)
  target_link_libraries(dynsyslab_colbrook_reference_tests PRIVATE dynsyslab_workflow dynsyslab_numerics)
  target_compile_definitions(dynsyslab_colbrook_reference_tests PRIVATE
    DYNSYSLAB_REFERENCE_DATA="${CMAKE_CURRENT_SOURCE_DIR}/tests/data")
  add_test(NAME dynsyslab_colbrook_reference_tests COMMAND dynsyslab_colbrook_reference_tests)
  set_tests_properties(dynsyslab_colbrook_reference_tests PROPERTIES
    LABELS "native;workflow;koopman;spectral;reference")
  add_executable(dynsyslab_rpe_hdmd_reference_tests tests/rpe_hdmd_reference_tests.cpp)
  target_link_libraries(dynsyslab_rpe_hdmd_reference_tests PRIVATE dynsyslab_workflow dynsyslab_numerics)
  add_test(NAME dynsyslab_rpe_hdmd_reference_tests COMMAND dynsyslab_rpe_hdmd_reference_tests)
  set_tests_properties(dynsyslab_rpe_hdmd_reference_tests PROPERTIES
    LABELS "native;workflow;rpe;hdmd;reference")
  add_executable(dynsyslab_linear_systems_postprocess_tests
    tests/linear_systems_postprocess_tests.cpp)
  target_link_libraries(dynsyslab_linear_systems_postprocess_tests PRIVATE dynsyslab_workflow)
  add_test(NAME dynsyslab_linear_systems_postprocess_tests
    COMMAND dynsyslab_linear_systems_postprocess_tests)
  set_tests_properties(dynsyslab_linear_systems_postprocess_tests PROPERTIES
    LABELS "native;workflow;postprocess;linear-systems")
  add_executable(dynsyslab_empirical_observability_postprocess_tests
    tests/empirical_observability_postprocess_tests.cpp)
  target_link_libraries(dynsyslab_empirical_observability_postprocess_tests PRIVATE dynsyslab_workflow)
  add_test(NAME dynsyslab_empirical_observability_postprocess_tests
    COMMAND dynsyslab_empirical_observability_postprocess_tests)
  set_tests_properties(dynsyslab_empirical_observability_postprocess_tests PROPERTIES
    LABELS "native;workflow;postprocess;linear-systems")
  add_executable(dynsyslab_bilinear_certificate_tests tests/bilinear_certificate_tests.cpp)
  target_link_libraries(dynsyslab_bilinear_certificate_tests PRIVATE dynsyslab_workflow dynsyslab_numerics)
  add_test(NAME dynsyslab_bilinear_certificate_tests COMMAND dynsyslab_bilinear_certificate_tests)
  set_tests_properties(dynsyslab_bilinear_certificate_tests PROPERTIES
    LABELS "native;workflow;postprocess;krom;control")
  add_executable(dynsyslab_bilinear_workflow_tests tests/bilinear_workflow_tests.cpp)
  target_link_libraries(dynsyslab_bilinear_workflow_tests PRIVATE dynsyslab_workflow dynsyslab_numerics)
  add_test(NAME dynsyslab_bilinear_workflow_tests COMMAND dynsyslab_bilinear_workflow_tests)
  set_tests_properties(dynsyslab_bilinear_workflow_tests PROPERTIES
    LABELS "native;workflow;postprocess;krom;control")
  add_executable(dynsyslab_mode_selection_tests tests/mode_selection_tests.cpp)
  target_link_libraries(dynsyslab_mode_selection_tests PRIVATE dynsyslab_workflow)
  add_test(NAME dynsyslab_mode_selection_tests COMMAND dynsyslab_mode_selection_tests)
  set_tests_properties(dynsyslab_mode_selection_tests PROPERTIES
    LABELS "native;workflow;postprocess;mode-selection")
  add_executable(dynsyslab_boundary_postprocess_tests
    tests/boundary_postprocess_tests.cpp)
  target_link_libraries(dynsyslab_boundary_postprocess_tests PRIVATE dynsyslab_workflow)
  add_test(NAME dynsyslab_boundary_postprocess_tests
    COMMAND dynsyslab_boundary_postprocess_tests)
  set_tests_properties(dynsyslab_boundary_postprocess_tests PROPERTIES
    LABELS "native;workflow;postprocess;krom;boundary")
  add_executable(dynsyslab_eigenvalue_optimization_postprocess_tests
    tests/eigenvalue_optimization_postprocess_tests.cpp)
  target_link_libraries(dynsyslab_eigenvalue_optimization_postprocess_tests
    PRIVATE dynsyslab_workflow)
  add_test(NAME dynsyslab_eigenvalue_optimization_postprocess_tests
    COMMAND dynsyslab_eigenvalue_optimization_postprocess_tests)
  set_tests_properties(dynsyslab_eigenvalue_optimization_postprocess_tests PROPERTIES
    LABELS "native;workflow;postprocess;koopman;spectrum")
  set(_dsl_case_study_root "${CMAKE_CURRENT_SOURCE_DIR}/../..")
  # Private parity test: runs only where the author's local K-ROM archive exists (never in a
  # public clone or a release), so the archive is located by pattern rather than named.
  file(GLOB _dsl_duffing_case_archive
    "${_dsl_case_study_root}/archives/derived/koopman-krom-extract/duffing-reference-free-training-30-exact-ics__*__f5b4dab7.np[z]")
  if(DYNSYSLAB_ENABLE_ARCHIVE_IO AND DYNSYSLAB_ENABLE_LINALG
      AND DYNSYSLAB_ENABLE_SUNDIALS AND EXISTS "${_dsl_duffing_case_archive}")
    add_executable(dynsyslab_klqr_case_study_parity_tests
      tests/klqr_case_study_parity_tests.cpp)
    target_link_libraries(dynsyslab_klqr_case_study_parity_tests
      PRIVATE dynsyslab_workflow)
    add_test(NAME dynsyslab_klqr_case_study_parity_tests
      COMMAND dynsyslab_klqr_case_study_parity_tests "${_dsl_case_study_root}")
    set_tests_properties(dynsyslab_klqr_case_study_parity_tests PROPERTIES
      LABELS "native;parity;klqr;case-study" TIMEOUT 120
      WORKING_DIRECTORY "${_dsl_case_study_root}")
  endif()
endif()

add_library(dynsyslab_workflow_c_api SHARED src/WorkflowCApi.cpp)
target_include_directories(dynsyslab_workflow_c_api PUBLIC include)
target_compile_features(dynsyslab_workflow_c_api PUBLIC cxx_std_17)
target_compile_definitions(dynsyslab_workflow_c_api PRIVATE
  DYNSYSLAB_WORKFLOW_C_API_BUILD=1)
target_link_libraries(dynsyslab_workflow_c_api PRIVATE dynsyslab_workflow
  PUBLIC dynsyslab_c_api)
if(MSVC)
  target_compile_options(dynsyslab_workflow_c_api PRIVATE /W4 /permissive- /EHsc /utf-8)
endif()

install(TARGETS dynsyslab_workflow dynsyslab_workflow_c_api dynsyslab_queue
  ARCHIVE DESTINATION lib LIBRARY DESTINATION bin
  RUNTIME DESTINATION bin)
install(DIRECTORY include/ DESTINATION include)

if(BUILD_TESTING AND DYNSYSLAB_BUILD_TESTS)
  add_executable(dynsyslab_krom_audit_query_tests tests/krom_audit_query_tests.cpp)
  target_link_libraries(dynsyslab_krom_audit_query_tests PRIVATE dynsyslab_workflow dynsyslab_workflow_c_api)
  add_test(NAME dynsyslab_krom_audit_query_tests COMMAND dynsyslab_krom_audit_query_tests)
endif()
