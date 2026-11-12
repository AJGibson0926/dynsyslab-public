# Native TimeSeries and Fourier contract

## Decision

`dynsyslab::core::TimeSeries` is the native C++ representation of one sampled
trajectory, control history, target history, or scalar observable. Python is
not required. The type owns immutable samples when data comes from an archive
and can expose an immutable zero-copy view when the solver already owns the
buffers.

Fourier analysis is an optional sidecar, not part of solution generation or
archive publication. This preserves fast solves, avoids redundant stored data,
and gives the CLI, Viewer, and future native algorithms one implementation.

## Ownership and layout

A `TimeSeries` contains:

- one increasing time vector of `N` samples;
- a row-major/interleaved `N x M` value matrix, where `M` is the component count;
- uniform-time-grid classification and, when valid, `deltaTime`;
- at most one replaceable in-memory Fourier cache.

Owning construction moves the source vectors into shared immutable storage.
`TimeSeries::view` allocates no sample-sized buffer; its source time/value
buffers must outlive the view. `NativeSolution::trajectorySeries`,
`NativeSolution::controlSeries`, and `TrajectoryBatch::trajectorySeries`
return these views directly over existing solver storage. Concrete matrices and
trajectory samples remain runtime objects and are never baked into a compile.

## Lazy Fourier sidecar

`series.fourier().transform(options)` is the only operation that performs a
transform. Construction, integration, catalog discovery, archive loading, and
ordinary time/phase display do not run it.

The implementation uses radix-2 Cooley--Tukey FFTs and Bluestein's algorithm
for arbitrary sample counts, preserving O(N log N) complexity without a new
runtime library. Options currently include:

- one-sided or two-sided output;
- optional mean removal;
- rectangular, Hann, Hamming, or Blackman windows.

The result exposes frequencies, complex coefficients, amplitude, power, phase,
and dominant-bin/frequency queries for every component. Window amplitudes are
corrected for coherent gain. Fourier input must use a uniform time grid.

Only the latest option set is cached. Repeating it returns the same immutable
result; requesting different options replaces the cache; `clear()` releases
it. No work or memory is retained unless a client explicitly requests the
analysis.

## Archive boundary

`storeTimeSeries` writes only authoritative time and value arrays.
`loadTimeSeries` materializes those arrays as an owning native object. Fourier
frequencies, coefficients, amplitudes, powers, and phases are deliberately
excluded from solution NPZ records because they are reproducible descriptions
of the stored samples.

Native solver archives declare `time_series_format` as
`dynsyslab.TimeSeries/v2`. Existing array names remain readable through the
archive adapter, so legacy data does not require rewriting.

## Viewer and CLI

Viewer's catalog contains compact metadata and archive pointers only. It does
not populate every archived trajectory at startup. Once the user selects a
frequency-domain view, Viewer loads only that record, applies its documented
coordinate transformation, and requests the native Fourier sidecar. Long
spectra are peak-preserving reduced for drawing only after the full FFT.

Available trajectory routes are:

- **Time** for ordinary sample histories;
- **Fourier** for component amplitude spectra;
- **Control frequency spectra** for saved controller histories, using the
  control coordinate-system metadata rather than state scaling;
- **Time + frequency** for both descriptions;
- overview/control routes, ordered control history, state time series, then
  phase space where those plots coexist.

The compiled CLI provides the same analysis as versioned JSON:

```powershell
dynsyslab-cli fourier trajectory.npz --component 0 --window hann
dynsyslab-cli fourier trajectory.npz --component 1 --keep-mean --two-sided
```

The command is read-only and requires no Python interpreter.

## Performance invariants

- No Fourier work occurs in a numerical RHS or solver output callback.
- Solver-to-TimeSeries access is zero-copy.
- Archive output adds no sample copy before NPZ encoding.
- Catalog startup reads no trajectory arrays.
- Explicit Fourier work is O(N log N), not a bounded O(N^2) DFT.
- Display decimation never changes the transform used for spectral analysis.
- Derived Fourier arrays never increase primary archive size.

Focused native tests cover arbitrary-length transforms, amplitude/frequency
accuracy, cache behavior, nonuniform-grid rejection, zero-copy views, archive
round trips, solver publication, Viewer route ordering, and catalog discovery.
