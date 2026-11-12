---
name: dynsys-delay-identification
description: Fit and audit H-SVD delay-coordinate dynamics with the standalone DynSysLab miscellaneous computes toolkit. Use for direct modal SINDy, structured tap models, rank or delay sweeps, and free forecast comparisons.
---

Use [the toolkit guide](../../../TOOLKIT.md) and the standalone `toolkit.py study` entry point. The [example](../../../examples/make_demo.py) generates a small study that can run without scientific archives.

- Translate the question into a delay geometry, retained rank, integration clock, and physical observable. More taps at fixed physical span and higher rank at fixed taps are different experiments. Capture the lag order, channel units, signed coordinates, and training scales.
- For a request to fit ξ′=f(ξ) from scratch, use `experiment: direct` and projected derivative labels. Declare instantaneous, weak-integral and trajectory-calibrated fits separately. A previous law can be a plotted baseline. Requests for inherited endpoint dynamics use the explicit structured recipe and saved model provenance.
- Split whole independent records before fitting the chart, dictionary, support, or readout. Use a common usable evaluation clock. Record any simulated prehistory and check a restart with entirely recorded history when startup could change the result.
- State the target of the regression explicitly. For an observation law returning present channel derivatives from a delay state, score it once on fresh measured windows and again in an autonomous delay rollout. Seed the complete initial history and build every subsequent input from the model's own past. The [side-car guide](../../../TOOLKIT.md#sensor-and-delay-reconstruction) documents the native solver and fractional lags. For waveform calibration freeze the STLSQ support, tune coefficients against training records, and keep held-out curves out of selection.
- Compare measured-window decoding, derivative fit on true states, and single-initial-window free rollouts. For value/rate input channels, inspect both the decoded rate and the derivative of decoded value. Use delay consistency and partial coverage when the present output looks convincing.
- Export the selected self-contained model. Use `toolkit.py replay` for a new initial condition and `equations` for the signed-state law and its provenance. Use `toolkit.py plot` with the full reference clock for reference-scale overlays and expanded model-scale pages. Keep the clock-aligned scores beside phase geometry.
