# dynsyslab archives

This is an independent data-and-archive package. It does not contain equations
of motion or numerical integration schemes.

All new numeric archives use lossless ZIP-LZMA compression inside the normal
NPZ container, so they remain directly readable with `numpy.load`.

- `solutions/<system>/` contains one atomic compressed numerical result for
  every completed solver call.
- `history/<system>/<run>/` contains optional streaming checkpoints and a
  manifest, preserved for diagnosis if a solve fails or is interrupted.
- `store.py` is the narrow archive-writing interface used by `solver/`.

Completed solution metadata includes `forcing_type` (`"forced"` or `"free"`)
and a compact `forcing` description. Built-in driver values are regenerated from
their parameters rather than duplicated sample-by-sample; only non-regenerable
external control values are saved. Hierarchies automatically deduplicate a
shared final time mesh into `time_grids/`, retaining each trajectory's raw
adaptive mesh only when it is necessary to preserve an interpolated solve.

Future indexing, sorting, inspection, and visualization tools should be added
to this package without coupling them to either `dynsyslib/` or `solver/`.
