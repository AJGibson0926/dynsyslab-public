# Trajectory resolution sidecar

This tool is available through the [miscellaneous computes side-car](../miscellaneous_computes/README.md) and can also be built independently.

Use this standalone tool to check whether apparent corners in a phase portrait come from sampling or drawing. It accepts any finite, real NPZ trajectory with a strictly increasing time vector and a sample-by-coordinate matrix. Choose two columns for a phase portrait. The clock comparison accepts any number of columns, including a scalar observable.

From the source root, build the C++ vertex selector and run its regression checks with a Python interpreter that has NumPy.

```powershell
$side = "tools/experiments/trajectory_resolution"
$build = "build/trajectory-resolution"
cmake -S $side -B $build -G Ninja -DCMAKE_BUILD_TYPE=Release "-DPython3_EXECUTABLE=$((Get-Command python).Source)"
cmake --build $build
ctest --test-dir $build --output-on-failure
$lib = "$build/trajectory_resolution.dll"
python "$side/export_vertices.py" --library $lib --source dense.npz --out display.npz --axes 0 1 --tolerance 0.00002
python "$side/resolution.py" dense.npz display.npz --sample-values state --window 0.76 0.78
python "$side/comparison_plot.py" --reference dense.npz --sample display.npz --sample-values state --window 0.76 0.78 --out comparison.png
```

The export stores both coordinate values, original row indices, timestamps, source path, sample count and tolerance. Use `--time-key` and `--state-key` for other archive schemas. The comparison commands accept separate keys for each archive. To audit an export of columns other than the default first two, pass a reference with those same two columns in the same order. The plot accepts custom axis labels, time label and legend text.

The selector keeps the endpoints and exact global coordinate extrema. It checks every intervening source sample against the straight segment between two retained samples, at the same fraction of the source-index interval. The largest allowed discrepancy in each coordinate is the requested fraction of that coordinate’s range. It subdivides failing segments. This bound describes the supplied discrete path. For irregular time grids, the separate clock audit evaluates time interpolation explicitly. `resolution.py` reports maximum absolute error, its time, sample counts and RMS on the overlapping clocks only. A finer independent solve can serve as a numerical reference.

The Python functions `compare_clocks`, `DisplayVertices`, and `observable_panels` can be used without the receiver study. The latter prepares separately scaled, labeled panels from arbitrary multichannel observables in base SI units.
