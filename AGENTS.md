# DynSysLab Engineering Policy

## First session with an agent

If the file `.dynsyslab/agent-welcomed` does not exist in this checkout, this is the user's
first agent session here: begin your first reply with a short greeting and the summary and
overview from `docs/AGENT_WELCOME.md` (in your own words, a few lines plus the "Get started"
steps), ask what they would like to do, and then create `.dynsyslab/agent-welcomed`
(git-ignored) so later sessions skip the welcome. Never repeat it once the file exists.
After installing a release for the user (`install.ps1`), also present that welcome and
overview in the chat. Start the desktop applications only when the user asks.

For installation or first-time use, read `docs/INSTALL.md` and use the verified
prebuilt release through `install.ps1`. No author-specific routing Hub or skill
is needed in a normal clone. Keep scientific archives, generated media, and
personal configuration local; commit source, small job specifications, reports,
and asset pointers. Do not mistake a local case-study archive path for bundled
sample data. `resources/examples/duffing-free-solve.json` is the installed
self-contained first-run example.

Public distribution is `AJGibson0926/dynsyslab-public`. The author's existing
`AJGibson0926/dynsyslab` repository retains private historical scientific media.
When publishing from that private checkout, export the reviewed tracked tree
into the public repository; never push or merge its historical branches into
the public repository. Ordinary changes made within a public clone can use
normal public branches. Keep the public installer pointed at the public release.

These rules apply to every new implementation and refactor in this repository.

Before editing, read `docs/START_HERE.md`, `ARCHITECTURE.md`, and the relevant
feature guide. Run `python -m dynsyslab doctor` to verify the checkout and
`python -m dynsyslab audit --format markdown` to inspect existing debt.

- Keep code lightweight, efficient, low-friction, flexible, and easy to remove or replace.
- Preserve clear architectural boundaries. UI packages compose public services; they do not duplicate numerical algorithms or own scientific data.
- Put reusable algorithms in the independently compilable C++ core and expose practical, documented C++ API and compiled CLI gates before coupling them to a UI. Python bindings and Python CLI helpers must remain optional thin clients over the native core.
- Prefer small modules with a single responsibility and explicit versioned data contracts at process boundaries.
- No source file may exceed 1,000 lines. Treat 500 lines as the normal design ceiling and split a module before crossing it unless a documented reason makes that harmful.
- Reuse the canonical archive, hierarchy, system, control, post-processing, Viewer catalog, and shared native UI components instead of creating parallel stores or frameworks.
- Preserve scientific provenance, coordinate-system metadata, archive immutability, and the Batch > System > Cluster > Trajectory hierarchy.
- Add focused tests and concise architecture/user documentation for new public behavior.
