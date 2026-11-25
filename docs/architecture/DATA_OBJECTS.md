# Data objects

Processing methods are black boxes: data objects go in, new data objects come out, and the method
keeps nothing. An archive holds one data object: a time series, eigenfunctions, reconstructions, a
Koopman model, a vector field, a field, a report. There is no "H-DMD archive" or "SINDy archive".
How a method made an object is a tag on the object, not where it lives.

## The envelope

Every archive a processing step writes carries `metadata_json["data_object"]`
(`core/DataObject.hpp`, schema `dynsyslab.data-object/v1`):

```json
{"schema": "dynsyslab.data-object/v1",
 "kind": "eigenfunctions",
 "produced_by": {"method": "hdmd", "settings": {...},
                 "inputs": [{"path": "...", "sha256": "...", "role": "input"}]},
 "properties": {}}
```

- `kind` names the contract of the archive's arrays.
- `produced_by` records the method, its runtime settings, and its inputs by path and SHA-256.
- Archives are filed by kind under `<archive root>/objects/<kind>/`. File names still carry the
  method, as a readable hint only.

Only `publish()` and the dispatcher stamp the envelope
(`workflow/native/src/ObjectPublishing.hpp`). The dispatcher records the running step's method and
settings. It also stamps any output that was written directly, once, before anything reads it.

## Kinds

| Kind | What it holds |
|---|---|
| `eigenfunctions` | A spectrum and its eigenfunctions, each stand-alone ([EIGENFUNCTION_OBJECTS.md](EIGENFUNCTION_OBJECTS.md)). |
| `reconstructions` | Observables written over one eigenfunctions object, g ≈ C ψ (+ offset). |
| `koopman_model` | Spectrum, eigenfunctions, readouts and input maps, ready for control or estimation. |
| `dynamical_model` | An identified vector field or rate law. |
| `time_series` | Trajectories and coordinate records, including recorded coordinate changes. |
| `coordinate_map` | Frozen bases and transforms (history bases, PCA). |
| `matrix` | A Hankel matrix a method read, when asked to keep it. |
| `field`, `figure` | Fields on domain grids, animations and figure datasets. |
| `report` | Analyses, validations, error tables and selections. |

## Composition

Objects refer to each other by path and SHA-256, never by copying.

- A `reconstructions` object names its basis in `metadata_json["reconstructions"]["basis"]`.
- `krom_extract` assembles a `koopman_model` from an eigenfunctions object, its source, and the
  optional setting `reconstructions` (a path or a list). `mergeReconstructions`
  (`core/KoopmanArchive.hpp`) checks the basis digest before attaching them.
- `loadKoopmanModel(path, binding)` attaches reconstructions listed in
  `ModelBindingOptions::reconstructions`.
- Viewer shows a reconstructions object on its basis. It shows an eigenfunctions object with the
  reconstructions written over it, found among the archive root's reconstructions objects.

## What H-DMD emits

H-DMD (Arbabi & Mezić, *SIAM J. Appl. Dyn. Syst.* 16, 2017; Brunton et al., *Nat. Commun.* 8, 2017)
publishes:

1. **eigenfunctions**, its primary output. It holds:
   - the spectrum;
   - the eigenfunction values at the training samples;
   - the shared interpolation (its domain, base points and interpolant);
   - each eigenfunction's own row of the window map.

   Per-eigenfunction properties: value error, Koopman residual, terminal modal amplitude and energy
   rank. The input's Hankel singular values are a property of the data the method read. The input
   series is referenced by its origin; no input values are copied.
2. **reconstructions**: the input channels (`reconstruction_matrix`), named readouts, and any
   catalog observables requested with `reconstruction_observables`, all written over those
   eigenfunctions.
3. **matrix**, only with `save_matrix`: the Hankel matrix and its shifted pair.

## What SINDy, DMD, EDMD and hand-given expansions emit

- **SINDy** (`sindy`, eigenfunction regression on a given spectrum; Brunton, Proctor & Kutz,
  *PNAS* 113, 2016; Kaiser, Kutz & Brunton, *Mach. Learn.: Sci. Technol.* 2, 2021) publishes:
  - **eigenfunctions**: the spectrum it was given and one polynomial expansion, each
    eigenfunction's coefficients its own row (`expansion_coefficients`). Per-eigenfunction
    properties are the training residuals. It stores no training states or coordinates; the
    series and the spectrum are referenced in its origin.
  - **reconstructions** of the requested observables, unless `"reconstruction": "none"`.

  `sindy_dynamics` publishes a `dynamical_model`.
- **DMD** (Schmid, *J. Fluid Mech.* 656, 2010; Tu et al., *J. Comput. Dyn.* 1, 2014) publishes
  **eigenfunctions** ψ = Φ⁺x as a linear expansion. The modes Φ and the reduced operator are
  properties.
- **EDMD** (Williams, Kevrekidis & Rowley, *J. Nonlinear Sci.* 25, 2015) publishes **eigenfunctions**
  in its own dictionary. The modes and the operator are properties.
- **`eigenfunction_expansion`** publishes **eigenfunctions** with continuous eigenvalues only.
  Hand-given rates have no sampling interval.

`loadKoopmanModel` builds a model from expansion eigenfunctions in any dictionary. With no
producing system, the model is coordinate-only.

## Pipelines

`${stage.archive}` and `${stage.output:METHOD}` name a method's primary object, its first.
`${stage.object:KIND}` names the stage's one object of that kind, for example
`"reconstructions": "${hdmd.object:reconstructions}"` in a `krom_extract` stage.

## Migration

- Archives written before data objects have no envelope. They are still read where their arrays
  are unchanged.
- Eigenfunction archives in the method layout are refused with a message to regenerate them. That
  layout means `domain__<name>`, or Hankel modes in place of `window_map`.
- Every processing step is now filed by kind and stamped. These still publish one combined
  `koopman_model` archive:
  - control models: DMDc and EDMDc;
  - K-ROM editing steps: `krom_extract`, recentring, readouts, input maps, control regression and
    mode-selection K-ROMs;
  - the spectral K-ROMs (mpEDMD, ResDMD and relatives);
  - boundary eigenfunctions;
  - the echo-history model;
  - coordinate SINDy on history or measured coordinates.

  The first two groups are models by nature. The spectral K-ROMs, boundary eigenfunctions and
  coordinate SINDy should still be split into eigenfunctions and reconstructions objects.
