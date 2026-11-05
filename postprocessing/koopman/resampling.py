"""Compatibility aliases for root-level postprocessing resampling."""
from ..resampling import ResamplingOptions as KoopmanResamplingOptions, resample as resample_for_koopman

__all__=["KoopmanResamplingOptions","resample_for_koopman"]
