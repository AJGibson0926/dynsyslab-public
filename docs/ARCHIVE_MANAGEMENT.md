# Archive management

Archives are NumPy `.npz` files (ZIP containers of `.npy` arrays plus `metadata_json`) under an
archive root. Completed archives are immutable: nothing here changes their scientific content.
Every command is part of the native CLI (`dynsyslab-cli`, or `dynsyslab archive|storage ...`
through the Python front door) and prints one JSON line.

## Where archives go

Without `--archive-root`, the CLI uses, in order:

1. the `DYNSYSLAB_ARCHIVE_DIR` environment variable;
2. the installation's `resources/runtime.json` (written by `install.ps1`);
3. `<checkout>/archives` when run from a source build;
4. `%LOCALAPPDATA%/DynSysLab/runtime.json`, then `%LOCALAPPDATA%/DynSysLab/archives`.

`dynsyslab-cli doctor` reports the choice and its source (`environment`, `install_runtime`, ...).
The desktop applications do not read `DYNSYSLAB_ARCHIVE_DIR` yet.

## Compression

Archives can be written with three lossless ZIP codecs, all readable by `numpy.load`:

| Codec | ZIP method | Use |
|---|---|---|
| `deflate` (default) | 8 | fast, compact |
| `lzma` | 14 (LZMA) | usually densest for smooth trajectories, slower to write |
| `stored` | 0 | no compression, fastest |

Choose the codec of every archive a command writes with `--compression lzma|deflate|stored`, or
for a whole session with `DYNSYSLAB_COMPRESSION`. Existing archives are re-encoded with

```powershell
dynsyslab-cli archive recompress --all-ages --compression lzma           # dry run: lists files
dynsyslab-cli archive recompress --all-ages --compression lzma --apply
dynsyslab-cli archive recompress results\a.npz results\b.npz --apply     # LZMA by default
```

Recompression writes the new file beside the old one, checks that every stored array and
`metadata_json` are byte-identical, and only then replaces the original. The archive index is
reconciled afterwards.

## Shared time grids

Many solutions of one campaign sample the same final time mesh. `share-time-grids` stores each mesh
used by two or more archives once, as `<root>/time_grids/<sha256>.npz`, and rewrites those
archives without their `time` array:

```powershell
dynsyslab-cli archive share-time-grids --section solutions --all-ages            # dry run
dynsyslab-cli archive share-time-grids --section solutions --all-ages --apply
```

The archive's metadata then names the grid (`time_grid_ref`, `time_grid_root`,
`time_grid_shared`, `time_grid_samples`). Every native reader restores `time` transparently; the
reference resolves against `time_grid_root`, then `archive_root`, then each folder above the archive,
so an archive tree that was copied elsewhere still reads. Archives written by the earlier Python
store with the same keys read the same way. Rewritten archives use the codec of the run
(`--compression`, `DYNSYSLAB_COMPRESSION`, else deflate). Raw adaptive meshes (`raw_time`) are evidence of a
particular solve and are never shared. Keep `time_grids/` with the archives it serves.

## Audit and deduplication

```powershell
dynsyslab-cli archive audit --all-ages --campaign <campaign-id> --report reports\
dynsyslab-cli archive deduplicate --all-ages --campaign <campaign-id>            # dry run
dynsyslab-cli archive deduplicate --all-ages --campaign <campaign-id> --apply
```

Scope options: `--days N` or `--all-ages`, `--section NAME` (repeatable, default `solutions`),
`--system NAME`, `--campaign ID` (archives whose metadata `campaign_id` matches),
`--include-derived`. `--report DIR` also writes `archive-audit__<UTC>.json` and `.md`.

The audit is read-only. It groups exact duplicates (identical `time` and `values`, compared byte
for byte) and repeated numerical plans. Deduplication moves only unreferenced exact duplicates into
archive-local trash with a recovery manifest. A duplicate that another archive or manifest refers
to is kept: rewriting the referring files would change immutable archives, so references are
never relinked.

## Campaign catalog

```powershell
dynsyslab-cli archive campaigns [--output catalog.json]
```

A read-only view of every campaign under the root, rebuilt from its manifests in `campaigns/` and
the `campaign_id` in solution metadata: label, kind, manifests, solution archives, and archives a
manifest lists that no longer exist (`dynsyslab.campaign-catalog/v1`). It stores nothing; the
queue's campaign registry is the stored record.

## Editing metadata

```powershell
dynsyslab-cli archive edit-metadata a.npz a-labelled.npz --metadata-json "{\"label\":\"run 3\"}"
```

Writes a new archive with the source's arrays and patched metadata; the source is unchanged. The
patch is a JSON object: each key is set, a `null` value removes it (JSON Merge Patch, RFC 7386, at
the top level). The edit is recorded in `metadata_edits` (source path and SHA-256, keys, UTC time)
and `metadata_revision` is incremented. Keys that define how an archive is read
(`recorded_archive`, `time_grid_*`, `run_provenance`, `metadata_edits`, `metadata_revision`)
cannot be patched.

## Storage locations, verified copy and relocation

```powershell
dynsyslab-cli storage register cold D:\dynsyslab-cold [--read-only] [--purpose solutions]
dynsyslab-cli storage list
dynsyslab-cli storage copy results\a.npz D:\backup [--remove-source]
dynsyslab-cli storage relocate <campaign-id> cold
dynsyslab-cli storage append <plan.json> <newly-completed.npz> ...
```

The registry (`management/storage_locations.json`, `dynsyslab.storage-registry/v1`) names folders
where archives may be kept; `primary` is the archive root. A verified copy hashes the copy with
SHA-256 against its source before publishing it, accepts an identical existing copy, never
overwrites a different file, and removes the source only on `--remove-source` after verification.

`storage relocate` creates a plan (`management/relocations/`, `dynsyslab.campaign-relocation/v1`)
for a registered writable location and copies every completed archive of the campaign (its
manifests, the hierarchy manifests and solutions they list, and solutions carrying its
`campaign_id`, and the shared time grids they reference), keeping their layout relative to the
archive root, so the copied tree reads on its own. Sources are kept. Files that
complete later are added with `storage append`; a file already written under the destination is
recorded without copying.

## References

- PKWARE, *APPNOTE.TXT – .ZIP File Format Specification*, version 6.3.10 (2022): compression
  methods 0 (stored), 8 (deflate) and 14 (LZMA).
- P. Deutsch, *DEFLATE Compressed Data Format Specification version 1.3*, RFC 1951 (1996).
- I. Pavlov, LZMA SDK (LZMA algorithm and format), https://www.7-zip.org/sdk.html.
- NIST, *Secure Hash Standard (SHS)*, FIPS PUB 180-4 (2015): SHA-256.
- P. Hoffman and J. Snell, *JSON Merge Patch*, RFC 7386 (2014).
- NumPy Enhancement Proposal 1, *A simple file format for NumPy arrays* (`.npy`/`.npz`).
