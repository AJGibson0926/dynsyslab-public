# Storage maintenance

This is host/workspace maintenance tooling, separate from the numerical pipeline.
It uses no scientific kernels. Inventory results, deletion manifests, archive
plans and paths identifying private data belong under ignored `research/`.

## Inventory

```powershell
python -u tools/maintenance/storage_inventory.py C:/ research/storage-review/c-before.json
```

Reads metadata only; skips directory reparse points; records access failures;
reports directory totals, large files, extensions and modification-age totals.
It never changes permissions or hydrates file contents. Logical bytes differ
from allocated space for sparse/compressed files. Modification age is not proof
of inactivity. Outputs must have fresh names.

## Reviewed cleanup

Prepare a private JSON policy with a `roots` array. Each rule has an absolute
`path`, `minAgeDays`, and a human-reviewed `reason`. Optional `extensions`
restricts selected files; optional `requireCMakeCache: true` checks a build root.
Use this only for files already established to be disposable. Cache/build names
alone are not evidence: build directories can contain irreplaceable results.

```powershell
./tools/maintenance/Invoke-StorageCleanup.ps1 -Mode Plan -PolicyPath research/storage-review/policy.json -PlanPath research/storage-review/plan.json
./tools/maintenance/Invoke-StorageCleanup.ps1 -Mode Apply -PlanPath research/storage-review/plan.json -ResultPath research/storage-review/result.json
```

Inspect the explicit file plan before applying. Application preflights containment
and allowed extensions, checks each file's size and modification time, and rejects
reparse ancestors. It deletes individual literal file paths, never recursively
deletes a directory. New, changed, locked and inaccessible files are preserved.
Source/configuration files and empty folders remain unless expressly selected.
The result records every deletion or skip. Live application profiles and databases
should be handled with their own maintenance/export mechanisms.

There is no filesystem transaction: do not apply concurrently with another task
that renames the same directories. A process interruption can stop before the
final JSON is written. Retain the plan and re-inventory in that case. Read free
disk space immediately before and after; background disk activity affects the
net change independently of the deletion manifest.

```powershell
python -m unittest discover -s tools/maintenance -p test_storage_cleanup.py -v
```

Tests cover unchanged planned-file deletion, preservation of fresh/results files,
changed files, whole-plan rejection on path escape, junction skipping, and a
directory replaced by a junction after planning.

## External archiving

An archive plan is distinct from a deletion policy. Record source paths, destination
relative paths, reason, size and preservation constraints. Identify the supplied
destination volume explicitly, copy first, verify hashes/counts/lengths, and only
then consider source retirement. Managed scientific archives require lineage and
reference checks through their owning module. Active SQLite, WSL and Docker files
need application-specific backup/export workflows; do not move them as ordinary
inactive files. Never follow junctions into unrelated trees during transfer.
