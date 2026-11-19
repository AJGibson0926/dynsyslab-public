# Observable Hankel matrices in the native pipeline

## Module boundary

| Module | Takes | Gives |
|---|---|---|
| Solver | system, initial condition, controller | state time series |
| Coordinate change (`coordinate_transform`) | any time series | the series in working coordinates: `identity`, an analytic `observable` map such as `echo_jet2`, `channels` (observables, their time derivatives, KROM observables) or `history` (Hankel-SVD coordinates, fitted here) |
| Hankel module | any time series and a window | the Hankel matrix: explicit (`buildHankelMatrix`) or matrix-free (`HankelOperator`), one row layout |
| Post-processing (H-DMD, Hankel SVD/POD, spectral DMD variants, SINDy) | a coordinate-change output | spectra, decompositions, eigenfunctions |

The Hankel module never evaluates physics. A channel is a column of the series
it receives, chosen by index or recorded column name. H-DMD and the history
basis both build their matrices through it; the history basis is fitted inside
the `history` coordinate change, never as a stage beside H-DMD. A
coordinate-change archive also stores the physical states at its samples
(`source_states`), so H-DMD keeps the physical state as its model state.

For older jobs and the workbench dialogs, a Hankel method that names physics
observables (as in the example below) is expanded by the orchestrator into an
explicit `coordinate_transform {"kind":"channels"}` step followed by the same
method reading that step's columns (`use_previous_output: true`). Consecutive
methods over the same channels share one coordinate change.

The `hankel_matrix` stage prepares one immutable matrix archive for each source trajectory. The `hankel_svd` (also `hankel_pod`) and `hdmd` (also `hankel_dmd`) stages can consume that exact archive. They can also construct their own inputs directly from the same settings. A job with multiple sources starts a fresh window sequence for each source.

## Run the RPE example

The [example job](../../examples/native/case-studies/rpe-observable-hankel.json) uses the archived RPE trajectory named in its `sources` list. Set this path to a local RPE solution if using another checkout. From the project root, run the compiled CLI.

```powershell
dynsyslab-cli plan examples/native/case-studies/rpe-observable-hankel.json
dynsyslab-cli run examples/native/case-studies/rpe-observable-hankel.json --archive-root archives
```

The job builds 512 windows with 32 consecutive samples per channel. Its three expanded channels are current echo, echo delayed by five analysis ticks, and the current time derivative of echo. Both analyses fit six modes from the same **96 × 512** matrix. To change retained rank on an existing matrix, supply the original trajectory as the source and set `hankel_archive` in the decomposition settings to the saved matrix path. The reader checks the source file's SHA-256. Within a job, `use_previous_hankel: true` selects the preceding matrix for the current source.

The workbench's H-DMD and Hankel POD settings accept comma-separated channel names or the JSON list below. Saved job specifications support the full common window contract.

```json
{
  "observables": [
    {"observable": "echo", "taps": [0, 5]},
    {"observable": "echo", "derivative_order": 1}
  ],
  "window": {"delays": 32, "columns": 512, "sample_stride": 1, "snapshot_stride": 1},
  "window_alignment": "end"
}
```

In this shorthand each entry selects a catalog name or alias, evaluated by the expanded coordinate change. An integer or recorded column name selects a column of the input series directly, including a measured rate supplied as its own column. Taps are nonnegative integer lags on the analysis grid. The default tap is zero. Channels expand in declaration order, then in the listed tap order. Optional `channel_scales` give one finite positive divisor per expanded channel. The default is one. Coincident channel declarations are rejected. Delaying an already included signal inside a dense window can produce repeated **rows** of the matrix; these retain their specified weighting and can reduce numerical rank.

## Indexing and alignment

Write the multichannel observation as $g(x_i)$ on the source clock with spacing Δt. Let $s$ be `sample_stride`, $a$ be `delay_stride`, and $b$ be `snapshot_stride`. All default to one. Let $L_c$ and $S_c$ be the tap lag and scale for expanded channel $c$, and $L=\max_c L_c$. With zero-based indices, $D$ delay blocks and $C$ expanded channels,

$$
 X_{dC+c,j}=g_c(x_{s(L+jb+da-L_c)})/S_c,
 \qquad
 Y_{dC+c,j}=g_c(x_{s(L+jb+da-L_c+1)})/S_c.
$$

Each column contains one window. Its chronological blocks run oldest to newest. For $a=b=1$, moving one block down gives exactly the same observation as moving one column right. Also $Y_{:,j}=X_{:,j+1}$ wherever both exist. Selecting every b-th column changes the training endpoints. The DMD prediction interval remains **one analysis tick**, $s\Delta t$. When $a\ne b$, the stored matrix is a sampled delay embedding without the adjacent-block anti-diagonal identity. The manifest records both strides. The extra taps change the measurements within each block and do not change the block advance.

The source ticks for column j's augmented start and end are $s(L+jb)$ and $s(L+jb+(D-1)a)$. An end-aligned snapshot is causal at its timestamp. The earliest raw sample in that column is $s j b$. The builder omits windows requiring unavailable history. Omitting `columns` or setting it to zero uses all complete windows. A paired matrix reserves one additional analysis tick for Y. A standalone SVD can use the final unpaired window too. Use a shared matrix stage when an identical set of columns matters. The matrix stage defaults to `shifted_pair: true`. Setting it to `false` builds X alone for SVD.

This path materializes dense matrices. X alone uses eight bytes per entry in memory, and a paired fit also needs Y and decomposition workspace. Set `window.columns` to bound the matrix size when working with long records. Embeddings too large to store, such as the pooled echo-history basis, use the matrix-free `HankelOperator` with the same layout: `coordinate_transform {"kind":"history"}` fits its left basis by randomized subspace iteration. The legacy `history_basis` stage uses the same operator.

The common `window` object uses explicit counts and defaults to end alignment. Legacy job fields keep their existing meanings: H-DMD's `matrix_size` counts pairs and `sampling_window` caps the delay count, while Hankel POD's `matrix_size` counts delay blocks. Legacy start-aligned jobs retain their timestamps. Choose either the common object or the legacy geometry fields.

## Observable derivatives

A system descriptor can declare a named first derivative and its domain of validity. Catalog lookups resolve aliases and can follow multiple links. For RPE the links are radius → radial velocity and echo → echo-dot → echo-double-dot. The named echo derivatives describe the zero-input drift on the model clock τ. Duffing links position → velocity and velocity → intrinsic force; the latter also describes zero-input dynamics.

`derivative_order` requests a derivative along the supplied trajectory. The resolver uses applicable catalog links, then analytic flow derivatives for a known drift or constant input, or an analytic observable Jacobian applied to available state velocities for a first derivative. Varying inputs can require input derivatives at higher orders. An unavailable derivative fails with an error. Requesting `"echo_dot"` by name explicitly selects the system's intrinsic observable, including when evaluated along a driven trajectory.

A system author registers a relationship using its canonical catalog names, for example `descriptor.observableTimeDerivatives["position"] = {"velocity", false};`. The boolean marks a link restricted to zero input. The catalog checks that both observables exist and have compatible dimensions and state coordinates. Aliases such as `r` and `echo` resolve through the same catalog.

To discover registered links without a state or a trajectory, pass the following query to `dynsyslab-cli observable query.json`.

```json
{"schema":"dynsyslab.observable-query/v1","operation":"list","system":{"name":"rpe"}}
```

## Mode and archive conventions

- **Hankel SVD / POD.** For $H=X=U\Sigma V^T$, `pod_modes` stores U (features × rank), `singular_values` stores Σ's diagonal, and `scores` stores $V\Sigma$ (windows × rank). Thus `pod_modes @ scores.T` reconstructs the retained approximation in scaled window coordinates. With `center: true`, add the saved per-feature `mean`. The historical arrays name the SVD of the transposed matrix: `left_singular_vectors` contains V and `right_singular_vectors` contains $U^T$. The sign of each singular-vector pair is arbitrary. SVD supplies measured temporal coordinates; it fits no evolution law.
- **Hankel DMD.** The existing projected DMD fits $\widetilde A=U_r^T Y V_r\Sigma_r^{-1}$. `modes` stores the projected DMD modes in scaled delay-feature space. `koopman_modes` and the named reconstruction arrays map the fitted coordinates to the aligned physical channels with channel scaling restored. The output contains discrete eigenvalues and their continuous-time logarithms on the recorded analysis clock. `eigenfunction_values` contains the fitted window coordinates at each X column plus the last Y column. These are estimates of Koopman coordinates. The archive's `training_time` gives their actual timestamps, including a shorter final interval when fitting columns were subsampled.
- **Provenance.** `metadata_json.hankel_window` records each resolved channel, derivative source, lag, scale, row order, strides, alignment, shape and analysis interval. `hankel_start_samples` and `hankel_end_samples` identify the original source rows. The matrix stage stores `hankel_matrix` and, for paired data, `shifted_hankel_matrix`. Direct-fit stages store these arrays only when requested with `save_matrix: true`.

The standalone matrix archive also carries `time` for every X column and `successor_time` for every stored Y column, both on the source trajectory clock.

## Code ownership

The system and observable catalog own derivative relationships. [ObservableDerivative.hpp](../../native/include/dynsyslab/core/ObservableDerivative.hpp) resolves their scientific meaning. [Hankel.hpp](../../native/include/dynsyslab/core/Hankel.hpp) builds a window or a shifted pair from sampled measurements without file or UI dependencies. Spectral DMD and POD consume that matrix through [Spectral.hpp](../../native/include/dynsyslab/core/Spectral.hpp) and [PCA.hpp](../../native/include/dynsyslab/core/PCA.hpp). The workflow evaluates channels, supplies the acquisition clock, and publishes matrices and decompositions through the archive API. The workbench edits the same JSON contract.
