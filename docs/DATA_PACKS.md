# Data packs

A data pack is a set of ordinary job files that regenerates a published set of results from
nothing. It holds parameters and instructions, never run data or images. Each pack is pinned to
one DynSysLab release, and runs offline with the native suite alone: no Python, no agent.

## Run a pack

Install the release the pack names, unzip the pack, then:

```powershell
$cli = "$env:LOCALAPPDATA\Programs\DynSysLab\bin\dynsyslab-cli.exe"
& $cli datapack validate .                          # files, hashes, portable job graph
& $cli datapack plan .                              # jobs, disk estimate
& $cli run jobs\00-reproduce-all.job.json --archive-root D:\pack-run
& $cli datapack verify . --archive-root D:\pack-run
```

Figures are written to `<archive-root>\job-groups\<pack>\figures`. Running the same command again
resumes. Each intermediate is deleted once its last consumer has finished and its digest is in the
ledger. The run refuses to start or continue below the pack's free-space floor (`--min-free-gb N`);
`--no-cleanup` keeps intermediates.

A different release refuses the pack and prints the right download link.
`--allow-version-mismatch` runs it anyway, labelled non-reproducing, and it can never verify.

## Verification

`verify` compares each job's canonical digests with the pack's sealed `expected.json`.

| Result | Meaning |
|---|---|
| `exact` | every array hash, or the PNG bytes, match |
| `numerically_equivalent` | hashes differ, but values agree within tolerance; never counted as verified |
| `failed` | a value, shape, array or image differs, or the job did not finish |

Array hashes cover type, shape and raw bytes. `metadata_json` is excluded because it holds times
and paths. Sealed digests keep each array's range, norm and eight samples, so a mismatch without
the original data can still show its size; the report marks that basis. With the original data,
`--reference-map MAP.json` compares every element and reports the worst one. PNGs are compared by
bytes, then by decoded pixels. The verdict is `verified_exact` only when every job is exact.

## Layout and authoring

```text
data-packs/<pack>/
  pack.json   pack.sha256   README.md
  jobs/00-reproduce-all.job.json     root job_group; names ../pack.json
  jobs/10-solves/ 20-koopman/ 30-input-maps/ 40-control/ 50-figure-data/ 60-figures/
```

`pack.json` (`dynsyslab.data-pack/v1`) holds the id, title, date, release pin, execution notes,
seeds, verification policy, figure list and the SHA-256 of every file. New schema versions get a
new name; older packs stay readable through `datapack inspect`. After editing a pack, run
`datapack rehash <pack>`. `datapack trace <archive.npz>` walks recorded lineage and lists gaps.

### Job kinds

- `job_group` runs child job files by relative path, in dependency order, with cycle and depth
  checks. A child names its archives (`"outputs":{"model":"archive"}`, or `"metadata:KEY"` for a
  path the stage recorded, such as `selected_krom`), and later jobs receive them
  through `"inputs":{"krom":"${fit.outputs.model}"}` and read `${inputs.krom}`. `${params.NAME}`
  passes parameters. A failed job blocks only its dependents. The ledger records every resolution.
- `viewer_figure` renders one PNG with the compiled Viewer CLI from complete JSON settings. A shared
  `style` scales fonts and line widths; `axis_limits.reference` takes limits from a declared render.

### Post-processing methods

- `krom_recenter`: re-expands polynomial eigenfunctions about `target_state` so psi(x*) = 0.
- `krom_readout`: adds a named readout; `lifted_split` is [Re psi, Im psi] with conjugate pairing.
- `reference_signal`: samples a target and its inverse-dynamics drive from the system's own RHS.
- `control_score`: tracking and drive metrics by window, with a declared selection rule.
- `figure_dataset`: selects, interpolates and scales channels into one plotting table.

The Hankel basis draws its start matrix from `mt19937_64+as241-inverse-cdf/v1`, using Wichura's
normal quantile (Algorithm AS 241, Applied Statistics 37(3), 477-484, 1988), so it no longer depends
on a standard library's `std::normal_distribution`.

## Release order

CI builds the release binary, so expected digests must come from that binary:

1. Commit, tag `vX.Y.Z` and push. CI builds and publishes the release.
2. Download and unzip that release's ZIP.
3. `tools/Build-DataPackAssets.ps1 -Cli <unzipped>\bin\dynsyslab-cli.exe -OutputDirectory <dir>`.
4. `gh release upload vX.Y.Z <dir>\DynSysLab-datapack-*.zip <dir>\SHA256SUMS-datapacks.txt`.
