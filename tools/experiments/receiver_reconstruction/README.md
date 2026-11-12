# Receiver reconstruction sidecars

These tools belong to the [miscellaneous computes side-car](../miscellaneous_computes/README.md). Its standalone build can also build the receiver and trajectory-resolution commands.

These standalone tools turn a saved spherical-bubble trajectory into synchronized receiver pressures, estimate state from labeled examples, and measure error by oscillation. They use the native RPE observable catalog and the standard NPZ format. Build this directory separately from the DynSysLab workflow and applications.

## Native commands

| Command | Inputs | Output |
| --- | --- | --- |
| `receiver_sidecar signals SOURCE.npz GEOMETRY.json OUT.npz` | One unforced RPE physical-state archive in native dimensionless coordinates, receiver positions in metres, sound speed in m/s, integer output stride | Receiver pressures in Pa, clock, reference state, delays, ring boundaries and source metadata |
| `receiver_sidecar estimate TRAIN.npz QUERY.npz K OUT.npz` | Training `features` (N×D), `state` (N×S), positive `scales` (D). Queries supply only `features` (M×D) and `time` (M) | State estimates (M×S), the exact K training row indices and weights per estimate, nearest distance, clock and scales |
| `receiver_sidecar rings REFERENCE.npz ESTIMATE.npz OUT.npz` | Matching time arrays, reference `state`, consecutive `ring_bounds`, optional `error_center` and estimated `prediction` | Per-ring counts, partial-window flags, error and reference sums of squares, componentwise relative RMS percentages |

The lookup library accepts other observation and state dimensions. It uses an exact nearest-neighbor search with deterministic row-index tie breaking. For each feature \(j\) it transforms the measurement to \(\operatorname{asinh}(y_j/s_j)\). For each query it finds the K nearest training rows in Euclidean distance and averages their **state labels** with weights proportional to \(1/\max(d,10^{-12})^2\). Training labels are used only by the inverse. Held-out labels are used only during scoring. Each estimate uses its own query. A preceding estimate or training timestamp does not constrain which rows are selected.

The `signals` clock is the nearest receiver's **emission time**. Column \(i\) contains pressure emitted \((D_i-D_\mathrm{near})/(cT_*)\) earlier, where \(T_*\) is seconds per unit of native time. The sidecar linearly interpolates each registered pressure term between its two bracketing source samples and starts once all paths have recorded history. Arrays are named after the pressure terms. For each array, `near_at_<term>` gives the near receiver's pressure at the **same emission** as the corresponding column. Configured positions can be in any order for this command. Sample states must remain inside the receiver range.

## Repeating the two-receiver comparison

`geometry.json` specifies the 5 mm and 6 mm receivers at right angles, 1480 m/s sound speed, and stride 40. `study.py` compares two distinct trajectories with matching parameters and geometry. Supply the nearer receiver first for this paired experiment. The scripts require NumPy and SciPy, with Matplotlib for figures. The native commands can run independently of Python.

From the source root, build against the matching native suite library. On Windows the DLL directory must be on PATH. The paths below correspond to this checkout's built native suite.

```powershell
$side = "tools/experiments/receiver_reconstruction"
$build = "build/receiver-reconstruction"
cmake -S $side -B $build -G Ninja -DCMAKE_BUILD_TYPE=Release "-DDYNSYSLAB_SUITE_BUILD=$((Resolve-Path build/native-suite-fast).Path)"
cmake --build $build
$env:PATH = (Resolve-Path build/global-jet/bin).Path + ";" + $env:PATH
ctest --test-dir $build --output-on-failure
$exe = (Resolve-Path "$build/receiver_sidecar.exe").Path
```

For each solver, supply separate physical training and held-out trajectories. The variables below denote local archive paths. Receiver geometry and fit settings must match between solver runs.

```powershell
& $exe signals $trainingArchive "$side/geometry.json" train-signals.npz
& $exe signals $heldOutArchive "$side/geometry.json" held-signals.npz
python "$side/study.py" --exe $exe --train train-signals.npz --held held-signals.npz --neighbors 8 --out rk4-study
# Repeat signals + study for the independent solver with fresh paths, producing bdf-study.
python "$side/diagnose.py" --exe $exe --study rk4-study --replicate bdf-study --output rebound.json
python "$side/figures.py" --study rk4-study --diagnosis rebound.json --out plots --library build/trajectory-resolution/trajectory_resolution.dll
```

Native outputs are created atomically at new paths. The study publishes a new directory after all comparisons finish. It saves the measurement definition, training and validation provenance, errors, digest of each scenario's inputs, a per-ring CSV, and diagnostics for replaying individual estimates. The figures command also writes a native Viewer catalog and payload with time series, physical \((r,r')\) and receiver-pressure portraits, correction overlays and neighbor weights.

Build the [trajectory resolution sidecar](../trajectory_resolution/README.md) for figures. Each physical reference curve reads the original solver archive. To also evaluate receiver pressure and the fitted inverse at **every** saved solver output in a selected interval, run the following command. Existing training examples, scaling and neighbor count are reused.

```powershell
python "$side/replay.py" --exe $exe --study rk4-study --window $start $stop --stride 1 --out early-native.npz
```

Repeat for the late window, then pass `--replays early-native.npz late-native.npz` to the figure command. Use the first eight intervals clipped to the common receiver clock and the last eight bounded intervals from the source ring boundaries. The plot verifies that both replays cover the exact native grids. The default whole-record estimates use the study clock. For a focused before-and-after comparison add `--detail-window START STOP --compare-vertices previous.npz`, where the previous export has `time` and two-column `state`. To show each acoustic term at each receiver, run `signals` with a geometry containing `"stride": 1` and `"window_tau": [START, STOP]` for the two source trajectories. Pass these archives to `figures.py --term-signals held-collapse.npz training-collapse.npz`. The window preserves earlier emission history for retarded evaluation.

## Comparing the observables

Every case starts with the monopole. The study adds the spherical term, dipole induction, or dipole radiation separately. It also runs the sequence spherical → induction → radiation. Only the **sum** at each receiver is given to the inverse. The same terms are included in training and evaluation. Identical feature matrices and scales across cases reuse an exact result, recorded under `computed_as`. This matters for the centered spherical model, whose two dipole fields vanish identically.

* **One pressure and two pressures.** The inputs are range-corrected totals \(D_i p_i\) taken at one receiver-clock time. They use a shared scale equal to the training RMS of the near channel. A nonzero range difference therefore includes echo history in the pressure pair.
* **Linear inverse.** The linearized free RPE around \((r,v)=(1,0)\) supplies the state matrix \(A\) and the intrinsic echo gradient \(H\). The rows \(H\exp(-A\delta_i)\) map the state deviations to the two range-corrected pressures. The inverse of that matrix estimates radius deviation and velocity. The spherical correction is quadratic in velocity, so its linearization at rest vanishes.
* **Ideal aligned history.** A separate three-feature experiment adds \(C=D_2p_2(t)-D_1p_1(t-\Delta)\) to the two simultaneous pressures, where \(\Delta=(D_2-D_1)/c\) in seconds. It uses only total pressures, including one additional near-receiver history sample. Under the spherical field model the common monopole cancels, leaving a multiple of the squared radial volume flux at the earlier emission. This feature uses its own training RMS scale. A contrast at roundoff level is omitted. The study also measures the error if the past near sample must be obtained by linear interpolation on the decimated receiver-output clock. A contrast-to-noise scale in Pa assumes independent equal-variance pressure noise and divides the RMS contrast by \(\sqrt{D_1^2+D_2^2}\). These are sensitivity diagnostics. The reconstruction inputs have zero added sensor noise.

## Rings, error and rebound attribution

A ring runs from one radius maximum to the next, identified by radial velocity changing from positive to nonpositive. An initial resting maximum counts as a boundary. Ring intervals are half-open. The first interval can be partially observed because of propagation delay, and the whole-record statistics include samples after the final complete ring. Relative RMS uses the RMS of \(r-1\) for radius and \(r'\) for velocity. Pool ring sums of squares before dividing when combining windows. Equal-weight samples correspond to a time average on a uniform clock.

`diagnose.py` selects the largest change in velocity estimation error after the fastest turn in the first rebound. It replays both adjoining estimates using their saved neighbor weights, groups those weights by the **training** ring, and compares with a second solver. A separate control keeps the held-out pressures fixed and halves the training sampling density. The resulting differences quantify how numerical repeatability and training-library sampling affect this particular inverse.
