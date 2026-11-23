# Spectral diagnostics tools

Native, stateless tools in `dynsyslab::math` that sit between fitted spectra, orbit tables and
sampled scalar functions. They read numbers and return reports. They own no archive and fit nothing.

| Module | Purpose |
| --- | --- |
| `EigenvalueLattice` | Flag eigenvalues near nonnegative-integer combinations of principal ones |
| `ShellPacket` | Predict and score rational-kernel spectral packets on a family of periodic orbits |
| `ReebGraph` | Critical points, persistence, merge tree, Reeb graph and separatrix levels on a 2D grid |
| `SpectralDiagnosticsQuery` | JSON adapter for the three modules |

## Eigenvalue lattice

Products of Koopman eigenfunctions are eigenfunctions, so eigenvalues of an autonomous system form
a semigroup (Mezic 2005; Budisic, Mohr & Mezic 2012). `eigenvalueLattice` tests
`lambda ~ sum_j m_j lambda_j` with `m_j >= 0`.

- **Coordinates.** Continuous eigenvalues are used directly. Discrete eigenvalues use `log(mu)`, so
  products become sums. A zero discrete eigenvalue is reported as `excluded`.
- **Aliasing.** Discrete angles are compared modulo `2*pi`. Continuous rates are compared modulo
  `2*pi/sample_interval` when `sample_interval` is given. `alias_turns` records the wrap.
- **Principal candidates.** Give `principal_indices`, or let a greedy pass pick them by ascending
  `|coordinate|`, skipping anything already explained, up to `maximum_principals`.
- **Matching.** The lowest order within tolerance wins, then the smallest detuning.
- **`chance_coverage`.** The fraction of the alias circle covered by all tolerance intervals. Near
  one, matches are expected by chance and carry little evidence.

`frequencyResonanceReport` applies this to Hankel or Spectral-KROM frequency lists. It folds signs
to `|omega|` and adds pairwise `m:n` candidates from the existing `resonanceCandidates`.

```json
{"schema":"dynsyslab.math-query/v1","operation":"eigenvalue_lattice",
 "clock":"discrete","eigenvalues":[[0.7,0.6],[-0.1,0.9]],
 "maximum_order":3,"maximum_principals":2,"tolerance":1e-3}
```

```json
{"schema":"dynsyslab.math-query/v1","operation":"frequency_resonance_report",
 "source":"hankel","frequencies":[0.7,1.3,2.0],"radians_per_sample":true,"tolerance":1e-6}
```

Results are **candidates, never rationality claims**. A match does not show that a mode is an
eigenfunction product, and the greedy generators need not be the true principal eigenvalues.

## Shell packets

For orbit coefficients `g_n(E)` and frequency `Omega(E)`, the predicted packet is

```text
psi(phi; E, theta) = sum_n w(n Omega dt; theta) g_n(E) exp(i n phi) exp(i n Omega c dt)
```

The weight is the rational kernel of Rigged DMD (Colbrook, Drysdale & Horning, arXiv:2405.00782):

```text
w(alpha;theta) = -(1/4 pi) sum_j [ r_j f(e^{i theta - i eps p_j})
                                 - conj(r_j) f(e^{i theta - i eps conj(p_j)}) ]
f(z) = (e^{i alpha} + z) / (e^{i alpha} - z)
p_j = i + (2j/(m+1) - 1),   sum_j r_j p_j^k = delta_{k0},  k = 0..m-1
```

It is the existing `RationalCircleKernel`, with `w(alpha;theta) = kernel(theta - alpha)`; the tests
check this against the literal formula. Order 1 is the Poisson kernel. Unitary spectral-measure
background follows mpEDMD (Colbrook, SIAM J. Numer. Anal. 61, 2023).

The module also provides:

- the smoothed density from `|g_n|^2`, through `smoothOrbitSpectra`;
- an alias catalogue solving `n Omega(E) dt = theta + 2 pi k` on a piecewise-linear `Omega(E)`;
- ridge width `h / (|n| dt |Omega'|)`, where `h` is the kernel's numerically computed half-width
  at half maximum;
- scores: relative L2 error, complex correlation, magnitude correlation and energy-profile
  correlation;
- a warning when a Krylov or delay target is one-sided rather than centred in its lag window;
- `n_last = c*w/h`. **This is an empirical heuristic.** The caller calibrates `c` and chooses the
  angular scale `w` (radians per sample, for example `Omega*dt`). It is not a resolution guarantee.

```json
{"schema":"dynsyslab.math-query/v1","operation":"shell_packet",
 "bandwidth":0.2,"order":2,"sample_interval":0.5,"angle":0.3,"target_delay":5,
 "maximum_harmonic":4,"delay_window":{"first_lag":0,"last_lag":10},
 "resolvable_constant":0.5,"angular_scale":0.6,
 "shells":[{"energy":0.0,"angular_frequency":1.0,"weight":0.5,"coefficients":[[1,1,0],[3,0.2,0]]},
           {"energy":1.0,"angular_frequency":0.8,"weight":0.5,"coefficients":[[1,1,0],[3,0.3,0]]}]}
```

```json
{"schema":"dynsyslab.math-query/v1","operation":"shell_score",
 "computed":[[1,0],[0,1]],"predicted":[[1,0],[0,0.9]],"shell_index":[0,1]}
```

Predictions cover only the supplied harmonics and table. Omitted harmonics, table interpolation
error and fitted-model error are not bounded.

## Reeb graph

`analyzeScalarGrid` works on the piecewise-linear interpolant over the Freudenthal triangulation of
a regular grid, with optional periodic axes. Ties break by vertex index.

- **Critical points** come from lower and upper link components (Banchoff 1970). Boundary vertices
  use their link path. Euler contributions sum to the domain's Euler characteristic.
- **Persistence** follows the lower-star filtration (Edelsbrunner, Letscher & Zomorodian 2002): H0 by
  union-find and the elder rule, H1/H2 by the same rule on the dual graph. Tests compare the pairs
  with a standard boundary-matrix reduction.
- **Merge tree** is the join tree of sublevel components (Carr, Snoeyink & Axen 2003).
- **Reeb graph** is built from components between consecutive critical levels. The join/split merge
  of Carr et al. holds only for simply connected domains, and periodic domains can have loops. Tests
  compare arc counts with directly counted contours at every regular level.
- **Separatrix levels** are interior saddles, plus boundary saddles where the graph truly branches.

`toFiniteGraph` returns a `FiniteGraph` for `graphTopology`, `graphDistances` and
`lowerStarPersistence`. In JSON, `edges` rows are `[from,to,length]`, as in the `graph` operation.

```json
{"schema":"dynsyslab.math-query/v1","operation":"reeb_graph",
 "nx":3,"ny":3,"values":[1,2,3,4,0,5,6,7,8],"periodic_x":false,"periodic_y":false,
 "minimum_persistence":0}
```

## Scope and limits

- The Reeb output describes the interpolant, not the underlying function. Under-resolved or noisy
  samples produce extra critical points, and the grid diagonal can bias near-degenerate saddles.
- `minimum_persistence` filters reported pairs only. The graphs are not simplified.
- Symmetric or flat data is resolved by index order, which can split a degenerate critical level
  into several nearby nodes.
- Grids are limited to 1,000,000 vertices. Reeb extraction rejects inputs over 20,000,000 slab pieces.

## References

- I. Mezic, Spectral properties of dynamical systems, model reduction and decompositions, Nonlinear
  Dynamics 41 (2005) 309-325.
- M. Budisic, R. Mohr, I. Mezic, Applied Koopmanism, Chaos 22 (2012) 047510.
- M. J. Colbrook, C. Drysdale, A. Horning, Rigged Dynamic Mode Decomposition, arXiv:2405.00782.
- M. J. Colbrook, The mpEDMD algorithm for data-driven computations of measure-preserving dynamical
  systems, SIAM J. Numer. Anal. 61 (2023).
- H. Edelsbrunner, D. Letscher, A. Zomorodian, Topological persistence and simplification, Discrete
  Comput. Geom. 28 (2002) 511-533.
- H. Carr, J. Snoeyink, U. Axen, Computing contour trees in all dimensions, Comput. Geom. 24 (2003)
  75-94.
- T. Banchoff, Critical points and curvature for embedded polyhedral surfaces, Amer. Math. Monthly
  77 (1970) 475-485.
