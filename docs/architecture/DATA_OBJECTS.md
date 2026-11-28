# Data objects

Processing methods are black boxes: data objects go in, new data objects come out, and the method
keeps nothing. An archive holds one data object: a time series, a Koopman reduced-order model
(krom), a vector field, a field, a report. There is no "H-DMD archive" or "SINDy archive".
How a method made an object is a tag on the object, not where it lives.

## The envelope

Every archive a processing step writes carries `metadata_json["data_object"]`
(`core/DataObject.hpp`, schema `dynsyslab.data-object/v1`):

```json
{"schema": "dynsyslab.data-object/v1",
 "kind": "krom",
 "produced_by": {"method": "hdmd", "settings": {...},
                 "inputs": [{"path": "...", "sha256": "...", "role": "input"}]},
 "properties": {}}
```

- `kind` names the contract of the archive's arrays.
- `produced_by` records the method, its runtime settings, and its inputs by path and SHA-256.
- Archives are filed by kind under `<archive root>/objects/<kind>/`. File names still carry the
  method, as a readable hint only.
- The method is the step's own name (`hdmd`, `sindy`, `krom_extract`, ...). The archive's
  top-level `method` is the same name, so the two never differ.

Every writer declares its object's kind (`core::object_kind`); nothing infers a kind from a name.
Archives are stamped when they are written: `publishKrom(root, source, model, name, metadata)`
writes and stamps a K-ROM, `publishObject(root, source, kind, name, archive)` files and stamps any
archive, `stampDataObject` stamps one a writer files itself, and `withDataObject` stamps metadata
handed to a core saver such as `saveKoopmanModel` (`workflow/native/src/ObjectPublishing.hpp`,
`PostprocessInternal.hpp`). Inputs default to the archives the metadata names as its sources. The
dispatcher only collects a step's objects, and refuses one without an envelope. Readers decide what
an archive is from its kind and contents, not from a method name.

## Kinds

| Kind | What it holds |
|---|---|
| `krom` | A Koopman reduced-order model: a spectrum whose point part are its eigenpairs (eigenvalue, eigenfunction with its samples and named representations, optional residual), its continuous part and optional pseudospectrum, optionally input maps B and readouts C. One that extends another references it as its parent. |
| `dynamical_model` | An identified vector field or rate law, including the DMDc and EDMDc controlled linear fits. |
| `time_series` | Trajectories and coordinate records, including recorded coordinate changes. |
| `coordinate_map` | Frozen bases and transforms (history bases, history transforms, PCA). |
| `matrix` | A Hankel matrix (from `hankel_matrix`, or from H-DMD with `save_matrix`); `hankel_archive` reads either. |
| `field`, `figure` | Fields on domain grids, animations and figure datasets. |
| `report` | Analyses, validations, error tables and selections. |

## The krom

A K-ROM **is** a spectrum and its eigenpairs. It may also carry an input map B (how inputs
drive the eigenfunctions) and readouts C (observables g ≈ C ψ + offset). Without B or C it is still
a model; they can be fitted later. In memory it is one class, `core::KoopmanModel`
(`core/KoopmanModel.hpp`); on disk it is one data object, schema `dynsyslab.krom/v1`
(`core/KoopmanArchive.hpp`).

- **Spectrum** (`core/SpectrumArchive.hpp`): the point part, `continuous_eigenvalues`,
  `discrete_eigenvalues`, `sample_interval`, `logarithm_branches`, `eigenvalue_domain`, and each
  eigenpair's residual when the method computed one, `eigenpair_residual` (NaN where a pair has
  none) with `eigenpair_residual_json` ({definition, scope} per pair); `mode_ids`, each pair's
  stable id when they are not 0..n-1; the pseudospectrum, `pseudospectrum_grid`,
  `pseudospectrum_minimum` and `pseudospectrum_json`; the continuous part, `spectral_bands_json`
  (bands; `continuous_spectrum_json` in older archives) and `spectral_measures_json` (each observable's
  measure on a band, as values). The functions on spectral sets (a band's generalized eigenfunctions,
  keys `band/<band>` and, for the level-set form, `band/<band>/level`; a measure's wave packets, key
  `measure/<band>/<observable>/<output>`; the pseudospectrum's minimizers, key `pseudospectrum`) are
  `spectral_function__<key>`, the
  observables' Koopman modes are the readouts (`reconstruction_*`), and a spectral
  fit's compressed operator is `compression_basis`, `compression_action` and
  `compression_readout__<name>`, both described in the krom block. See
  [EIGENFUNCTION_OBJECTS.md](EIGENFUNCTION_OBJECTS.md#eigenpairs-spectrum-and-pseudospectrum).
- **Samples**: `eigenfunction_values` (samples x modes), `training_time` and
  `training_sample_indices` (rows of the input series the origin points to). A fitted expansion
  may have none.
- **Representations**, each one multi-output function on its own domain, named and ordered; the
  first is the one the model evaluates (`KoopmanModel::useRepresentation`):
  `expansion_coefficients[__<name>]` (modes x terms; polynomial, DMD's linear map, an EDMD or
  spectral dictionary), `interpolation_points` (samples x domain dimension; its values are the
  samples), `window_map` (modes x window rows: the time-delay map of H-DMD). Descriptions are in
  `metadata_json["krom"]["representations"]`, with `order` and `default_representation`.
- **Readouts C**: `reconstruction_names` (newline-joined), `reconstruction__<key>` (outputs x modes),
  `reconstruction_offset__<key>`, `reconstruction_error__<key>` (rms, max, relative).
- **Input maps B**: `control_<name>_*`, `discrete_input_map_<i>` and `input_map_manifest_json`; each
  records the identity of the eigenfunctions it acts on (`koopmanModelIdentity`).
- **Diagnostics and provenance**: `spectral_diagnostics_json` (invariance, stochastic residuals and
  RKHS adjoints of a spectral fit; wave packets belong to their measure in the continuous spectrum), `dataset_references`, `binding_evidence`,
  `metadata_json["krom"]["origin"]` (method, settings, inputs) and `input_channels`.

Arrays a method adds (modes, operators, Hankel singular values) stay beside these as the producing
method's properties. Eigenpair residuals are not among them: they live on the spectrum, once.

### Parents

Objects refer to each other by path and SHA-256, never by copying. A step that adds to a krom it
read publishes a new krom whose `metadata_json["krom"]["parent"] = {path, sha256, modes}` names the
parent and the parent's modes it keeps, in order. The child stores only what it adds: a readout, an
input map, another representation, the residuals of its pairs (SINDy's generator residuals), or
nothing at all (a mode selection). Eigenvalues and samples always come from the root; the kept
modes keep their residuals.

`loadKoopmanModel(path)` reads any krom: it loads the parent (refused when it is missing or its
hash changed), keeps the listed modes, then lays the child's representations, readouts and maps
over them (one of the same name replaces the parent's). A loaded model remembers the krom it was
loaded from (`KoopmanModel::parent`); everything it has then is inherited, and what a step adds
afterwards is its own, so writing it produces the child automatically. `clearParent` writes a
whole krom instead (mode selection's conjugate-consistent rewrite does).

A krom without B gives a clear error when control asks for one ("this K-ROM has no input map B
named ..."). A time-delay K-ROM (its first representation a window map) evaluates on windows only;
asking it for a state value or a state-Jacobian map is refused.

| Step | Publishes |
|---|---|
| `hdmd`, `dmd`, `edmd`, history EDMD (`echo_history_koopman`), `koopman_spectral` (`resdmd`, `mpedmd`, `hermitian_dmd`, `tlsdmd`, `multdmd`, `stochastic_resdmd`, `specrkhs`), `boundary_eigenfunctions`, `eigenfunction_expansion`, coordinate SINDy, `sindy` from data | a whole krom |
| `sindy` given a krom, `eigenfunction_polynomial_degree` | a child adding its polynomial (the evaluated representation by default) |
| `observable_reconstruction` | a child adding a readout C |
| `input_map`, `control_regression`, `history_control_regression` | a child adding input maps B |
| `krom_extract`, `krom_readout`, K-ROM validation links | a child with the selected modes (their eigenpairs keep their residuals; and any refitted readouts) |
| `mode_selection` | a whole krom (its conjugate pairs are rewritten; a rewritten lower mode takes its partner's residual), with readouts and B |
| DMDc/EDMDc `export_krom` | a whole krom with its samples, B and consistency |
| `krom_recenter` | a whole krom (re-expanded eigenfunctions) |

## What H-DMD emits

H-DMD (Arbabi & Mezić, *SIAM J. Appl. Dyn. Syst.* 16, 2017; Brunton et al., *Nat. Commun.* 8, 2017)
publishes:

1. **krom**, its primary output. It holds:
   - the spectrum;
   - the eigenfunction values at the training samples;
   - the shared interpolation (its domain, base points and interpolant) and each eigenfunction's
     own row of the window map;
   - readouts C of the input channels (one per channel, and `reconstruction_matrix`) and of any
     catalog observables requested with `reconstruction_observables`.

   Its eigenpairs have no residual (H-DMD computes none). Per-eigenfunction properties: value
   error and Koopman residual of each representation, terminal modal amplitude and energy rank. The input's Hankel singular values are a property of the data the method read. The input
   series is referenced by its origin; no input values are copied.
2. **matrix**, only with `save_matrix`: the Hankel matrix and its shifted pair.

## What SINDy, DMD, EDMD, spectral methods and hand-given expansions emit

- **SINDy** (`sindy`, eigenfunction regression on a given spectrum; Brunton, Proctor & Kutz,
  *PNAS* 113, 2016; Kaiser, Kutz & Brunton, *Mach. Learn.: Sci. Technol.* 2, 2021) publishes a
  **krom** with one polynomial expansion (`expansion_coefficients`) and the readouts of the
  requested observables, unless `"reconstruction": "none"`. Given a krom (`krom`, `eigenfunctions`,
  `spectrum_archive` or the preceding step's object) it is that krom's child: the same
  eigenfunctions, one more representation, calibrated to their sample values. Each eigenpair
  carries its generator residual (`generator`, RMS over the training rows of the published,
  calibrated eigenfunction; with `make_default: false` the parent's residuals stay); the largest
  pointwise residual is the property `training_residual_maximum`. It stores no training states or
  coordinates.
  `sindy_dynamics` publishes a `dynamical_model`.
- **DMD** (Schmid, *J. Fluid Mech.* 656, 2010; Tu et al., *J. Comput. Dyn.* 1, 2014) publishes a
  **krom** with ψ = Φ⁺x as a linear expansion. The modes Φ and the reduced operator are properties.
- **EDMD** (Williams, Kevrekidis & Rowley, *J. Nonlinear Sci.* 25, 2015) publishes a **krom** in its
  own dictionary. The modes and the operator are properties.
- **Spectral K-ROMs** (`koopman_spectral`; ResDMD: Colbrook & Townsend, *Comm. Pure Appl. Math.* 77,
  2024) publish a **krom** with the frozen-dictionary expansion and an interpolation on the physical
  states at the training columns; the one that is a function of the state is evaluated first. Each
  eigenpair carries its one-step residual (`one-step`, `empirical-training`) of the frozen-dictionary
  expansion; a `residual_grid` becomes the spectrum's pseudospectrum; measures (mpEDMD, Rigged DMD,
  ResDMD; Colbrook, *SIAM J. Numer. Anal.* 61, 2023; Colbrook, Drysdale & Horning,
  arXiv:2405.00782) its spectral bands; the rest the krom's spectral diagnostics.
- **`eigenfunction_expansion`** publishes a **krom** with continuous eigenvalues only. Hand-given
  rates have no sampling interval.

A DMD or EDMD fit whose eigenfunctions cannot be formed as an expansion publishes its spectrum,
modes and operator as a `report`. With no producing system, a krom is coordinate-only.

## Viewer

The Viewer shows a child krom composed: the whole K-ROM it describes, with its ancestors' method
arrays (training samples, diagnostics) where it has none (`viewer/native/src/KromTraining.hpp`,
`composedKrom`). A krom whose evaluated representation is a function of the state opens in the
K-ROM dashboard.

## Experiments

`${stage.archive}` and `${stage.output:METHOD}` name a method's primary object, its first.
`${stage.object:KIND}` names the stage's one object of that kind, for example
`"${hdmd.object:krom}"`.

## Migration

- Archives written before data objects have no envelope. They are still read where their arrays
  are unchanged.
- The kinds `eigenfunctions`, `reconstructions` and `koopman_model` (schemas
  `dynsyslab.eigenfunction-objects/v1` and `dynsyslab.koopman-model/v3`) are gone; their archives
  are refused with a message to regenerate them. The step setting `reconstructions` (readouts
  objects to attach) is refused: readouts belong to the krom.
