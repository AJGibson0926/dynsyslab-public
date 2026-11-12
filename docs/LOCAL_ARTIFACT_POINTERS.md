# Local scientific and media artifacts

DynSysLab keeps generated scientific results and rendered media outside Git.
The repository versions source code, reproducible job configurations, narrative
reports, and lightweight provenance. It does not version new solver archives,
Viewer payloads, screenshots, plots, PDFs, or generated archive indexes.

## Local locations

- `archives/` is the canonical local scientific archive. Its solution,
  derived-data, campaign, management, and system-metadata products remain on
  the workstation and are addressed through archive identifiers recorded in
  reports and job configurations.
- `viewer_outputs/` contains local Viewer renders, comparison payloads,
  temporary Viewer archives, plots, and PDFs. Case-study documents refer to
  these products by relative path when an on-machine rendering is useful.
- `docs/reports/*/assets/` contains screenshots and generated diagnostics used
  by local HTML reports. The report source is versioned; the rendered assets
  remain local.
- `viewer_compute_archive/`, `test_reports/`, `logs/`, `manifests/`, and other
  runtime directories remain governed by the existing ignore policy.

These exclusions do not delete or relocate local data. A future publication
workflow may promote selected artifacts to an external object store or a
release attachment with checksums and stable URLs. Until then, the repository
contains only the pointers and reproduction inputs needed to locate or
regenerate them on the originating workstation.
