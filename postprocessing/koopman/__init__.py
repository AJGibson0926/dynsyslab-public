"""Coherent Koopman/DMD analysis API for archived trajectories."""
from .dmd import HDMD_DEFAULT_MATRIX_SIZE, HDMD_DEFAULT_SAMPLING_WINDOW, HDMD_DEFAULT_SVD_TOLERANCE, dmd, hankel_dmd, hankel_dmd_ensemble
from .edmd import edmd
from .measures import empirical_spectral_weight
from .dictionary import ObservableDictionary, PolynomialDictionary, RBFDictionary
from .result import KoopmanResult
from .pipeline import dmd_archive, dmd_processor, edmd_archive, hankel_dmd_archive, hankel_dmd_processor
from .validation import require_uniform_time_steps
from .validation import ResamplingDecisionRequired
from .resampling import KoopmanResamplingOptions, resample_for_koopman
from .spectrum import (DynamicalSystemReference, KoopmanEigenvalue, KoopmanEigenfunction,
                       KoopmanSpectrum, ScatteredEigenfunctionInterpolator)
from .koopman_sindy import KoopmanSINDyResult, koopman_with_sindy
from .error_analysis import (EigenfunctionErrorEvaluation, ErrorStatistics,
                             FreeTrajectoryDataset, SemigroupHorizonError,
                             evaluate_eigenfunction_errors)
from .eigenvalue_optimization import EigenvalueOptimizationResult, optimize_trajectory_eigenvalues
from .boundary_eigenfunctions import BoundaryEigenfunctionResult, fit_boundary_eigenfunctions
from .composable import (ComposableKoopmanResult, EigenfunctionStageResult,
                         EigenvalueStageResult, KoopmanPipelineSpecification,
                         construct_eigenfunctions, eigenvalues_from_hdmd,
                         provided_eigenvalues, run_composable_koopman_pipeline)
from .pipeline import koopman_with_sindy_archive
from .workflows import (DuffingControlledTrainingSOP, DuffingKoopmanSOP,
                        RPEControlledTrainingSOP, RPESymmetricPulseSOP, duffing_hdmd_sindy,
                        duffing_reference_initial_conditions,
                        generate_duffing_controlled_training,
                        generate_duffing_free_training,
                        generate_rpe_controlled_training,
                        generate_rpe_symmetric_pulse_training,
                        rpe_controlled_initial_conditions,
                        rpe_symmetric_pulse_anchors,
                        sample_duffing_controlled_inputs,
                        sample_rpe_controlled_inputs)
from .observable_spectra import (
    HDMDObservableComparison,
    ObservableSpectrumEntry,
    compare_hankel_observables,
    compare_hankel_observables_archive,
    load_observable_spectrum_comparison,
    observable_spectrum_comparison_report,
    save_observable_spectrum_comparison,
)

__all__=["dmd","hankel_dmd","hankel_dmd_ensemble","HDMD_DEFAULT_MATRIX_SIZE","HDMD_DEFAULT_SAMPLING_WINDOW","HDMD_DEFAULT_SVD_TOLERANCE","edmd","koopman_with_sindy","empirical_spectral_weight","ObservableDictionary","PolynomialDictionary","RBFDictionary","KoopmanResult","KoopmanSINDyResult","DynamicalSystemReference","KoopmanEigenvalue","KoopmanEigenfunction","KoopmanSpectrum","ScatteredEigenfunctionInterpolator","EigenfunctionErrorEvaluation","ErrorStatistics","FreeTrajectoryDataset","SemigroupHorizonError","evaluate_eigenfunction_errors","EigenvalueOptimizationResult","optimize_trajectory_eigenvalues","BoundaryEigenfunctionResult","fit_boundary_eigenfunctions","ComposableKoopmanResult","EigenfunctionStageResult","EigenvalueStageResult","KoopmanPipelineSpecification","construct_eigenfunctions","eigenvalues_from_hdmd","provided_eigenvalues","run_composable_koopman_pipeline","dmd_archive","hankel_dmd_archive","edmd_archive","koopman_with_sindy_archive","dmd_processor","hankel_dmd_processor","require_uniform_time_steps","ResamplingDecisionRequired","KoopmanResamplingOptions","resample_for_koopman","HDMDObservableComparison","ObservableSpectrumEntry","compare_hankel_observables","compare_hankel_observables_archive","load_observable_spectrum_comparison","observable_spectrum_comparison_report","save_observable_spectrum_comparison","DuffingControlledTrainingSOP","DuffingKoopmanSOP","RPEControlledTrainingSOP","RPESymmetricPulseSOP","duffing_hdmd_sindy","duffing_reference_initial_conditions","generate_duffing_controlled_training","generate_duffing_free_training","generate_rpe_controlled_training","generate_rpe_symmetric_pulse_training","rpe_controlled_initial_conditions","rpe_symmetric_pulse_anchors","sample_duffing_controlled_inputs","sample_rpe_controlled_inputs"]
