---
name: dynsys-sensor-reconstruction
description: Compare receiver configurations, delayed measurements, sensor corrections, and state readouts on common trajectories using the standalone DynSysLab miscellaneous computes side-car.
---

Use the [reconstruction manifest](../../../TOOLKIT.md#sensor-and-delay-reconstruction). For spherical-bubble pressure generation and correction terms, use the [existing receiver companion](../../../../receiver_reconstruction/README.md). Feed its physical channels into the generic comparison, preserving their units and emission-time convention.

Construct scenarios explicitly from channel, physical lag and scale. Use common complete history across scenarios. A change in receiver distance can change amplitude and arrival time; check how alignment and source-term normalization affect the comparison. A derivative estimated from a pair of pressure samples refers to the midpoint and has half-spacing latency.

Fit readouts on the labeled training trajectories. Compare affine and neighbor readouts on fresh held-out measurements. If autonomous forecasting is requested, initialize once and advance a learned map; report its advance interval and full or partial coverage. A raw-delay discrete map and an H-SVD rank-truncated discrete map use different input coordinates. Use a frozen modal chart for the latter.

Use matched per-ring metrics for transients. In Viewer, keep the reference, measurement-driven reconstruction and free forecast visibly identifiable and expose the collapse on its own time scale when needed.
