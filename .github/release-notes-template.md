### Install and run (Windows x64, no admin, no Python)

```powershell
Invoke-WebRequest https://github.com/AJGibson0926/dynsyslab-public/releases/download/{{TAG}}/install.ps1 -OutFile install.ps1
powershell -ExecutionPolicy Bypass -File .\install.ps1 -Version {{TAG}}
& "$env:LOCALAPPDATA\Programs\DynSysLab\bin\dynsyslab-cli.exe" quickstart
```

**Running DynSysLab means `dynsyslab-cli quickstart`.** It prints a welcome and overview of
DynSysLab, runs the bundled examples (a Duffing double well and the chaotic Lorenz system),
and reports their results as Markdown: a table of each run (samples, solver steps, dominant
frequency, final state), the archive paths, and PNG figures rendered by the Viewer CLI. The
installer runs it once at the end (skip with `-NoQuickstart`). Results are saved in
`Documents\DynSysLab\archives`. No window opens. Nothing is shipped pre-computed: every example
generates its own data.

**Chat with DynSysLab.** The installer also connects Claude Code and Codex (the `dynsyslab` MCP
server, a Claude Code start-up hook, and a workspace in `Documents\DynSysLab`; skip with
`-NoAgentSetup`, undo with `dynsyslab-cli agent-setup --remove`). From the next chat session on,
the agent opens with the DynSysLab Hub: a welcome and overview on your first login after an
install, then "What would you like to do?" with options. `dynsyslab-cli hub` shows the same.

The desktop applications (Solver/Post-processor and Viewer) are optional; open them from the
Start menu when you want them. More: [INSTALL.md](https://github.com/AJGibson0926/dynsyslab-public/blob/{{TAG}}/docs/INSTALL.md)
and the [overview](https://github.com/AJGibson0926/dynsyslab-public/blob/{{TAG}}/docs/AGENT_WELCOME.md).

Python client (optional, any OS): `python -m pip install git+https://github.com/AJGibson0926/dynsyslab-public@{{TAG}}`

### Data packs

Data packs are version-pinned job files that regenerate published data and figures from scratch,
offline, with this release alone. They hold parameters and instructions, never run data. Packs for
this release are attached as `DynSysLab-datapack-<name>-<release-date>-v<version>.zip`, with
checksums in `SHA256SUMS-datapacks.txt`. Unzip one and follow its `README.md`; see
[DATA_PACKS.md](https://github.com/AJGibson0926/dynsyslab-public/blob/{{TAG}}/docs/DATA_PACKS.md).
