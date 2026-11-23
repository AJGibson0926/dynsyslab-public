# RPE echo-only KLQE shadow start-up figures

Pack `rpe-echo-klqe-shadow-startup-2026-11-23`, pinned to DynSysLab v0.8.0. It holds instructions only: every result is generated on your machine.

Offline, with the installed v0.8.0 release:

```powershell
$cli = "$env:LOCALAPPDATA\Programs\DynSysLab\bin\dynsyslab-cli.exe"
& $cli datapack validate .
& $cli datapack plan .
& $cli run jobs\00-reproduce-all.job.json --archive-root D:\pack-run
& $cli datapack verify . --archive-root D:\pack-run
```

The run keeps at least 6 GB free (`--min-free-gb N`), deletes each intermediate after its last
consumer, and resumes if it is run again. Figures are written to
`<archive-root>/job-groups/rpe-echo-klqe-shadow-startup-2026-11-23/figures/`. See `docs/DATA_PACKS.md`.
