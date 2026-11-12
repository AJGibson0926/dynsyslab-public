---
name: dynsys-closure-witness
description: Investigate apparent phase-space crossings, derivative-jet closure and resolution sensitivity with DynSysLab. Use for branch witnesses, higher-coordinate lifts and routing a global symbolic closure question.
---

For sampled trajectories, use [crossings and refinement](../../../TOOLKIT.md#crossings-and-closure). For symbolic queries the same guide shows `toolkit.py closure` with a compatible native CLI. Search for numerical witnesses on the full scientific record. Record axis units/scales, both original segment indices, physical visit times and the separation tolerance.

Refine the two visits on independent high-resolution data or a continuous observation from the solver. Compare the scaled mismatch and the two velocities. Lift the visits into any additional coordinates before attributing a crossing of two plotted coordinates to the full proposed state. Match visits across solver tolerances or time grids with `witness_checks.match_visits`; retain unmatched candidates. Display simplification can create or remove apparent intersections.

A numerical witness can falsify a proposed single-valued derivative if the same state has incompatible derivatives, subject to the evidence for equality. Treat finite-precision equality and tests over a finite record as empirical. For a global question about analytic equations, use the existing [native closure guide](../../../../../../docs/architecture/COORDINATE_DYNAMICS_SINDY.md) and its v2 query for supported certificate classes. Report the scope, assumptions, witness or certificate, and unresolved orders explicitly.
