# Welcome to DynSysLab

*Shown by AI coding agents (Claude Code, Codex, Gemini CLI and others) the first time they
open this repository. Agents: present the greeting and overview below in your own words,
briefly, then ask what the user wants to do.*

---

**Hello, and welcome to DynSysLab.** It is a native C++ suite for simulating, analysing,
modelling and controlling dynamical systems: from equations, through data-driven Koopman
models, to estimators and controllers, with every result archived and every figure traceable
to its data.

## What it does

1. **Simulate** a library of ODE, delay, lattice and PDE systems (bubble dynamics/RPE, Duffing,
   Lorenz, pendula, Hodgkin–Huxley, Kuramoto, …), free or forced, single runs or campaigns.
2. **Analyse** trajectories: SVD/POD, Fourier spectra, phase space, Lyapunov screens.
3. **Identify models**: DMD family, Hankel/H-DMD, SINDy and spectral Koopman methods produce
   Koopman reduced-order models (K-ROMs) with eigenvalues, eigenfunctions, readouts and input maps.
4. **Estimate and control**: Koopman LQE (KLQE) and Koopman LQR/MPC, classical LQR, PID,
   integral action, and a linear-systems toolkit (gramians, balanced truncation, transfer
   functions, LQG, loop margins, observability measures).
5. **View** everything in the Viewer or headlessly with its CLI: time series, phase portraits,
   spectra, Bode/Nyquist, stem plots, multi-run comparisons.

Every computation is a versioned JSON job (`dynsyslab.compute-job/v1`) that can be validated,
planned, run and re-run; archives are immutable and record their provenance.

## Get started

- **Install (Windows x64, no compiler or Python needed):** run [`install.ps1`](../install.ps1),
  which downloads and checksum-verifies the latest
  [release](https://github.com/AJGibson0926/dynsyslab-public/releases/latest). Details:
  [INSTALL.md](INSTALL.md).
- **Python client (optional, any OS):**
  `pip install git+https://github.com/AJGibson0926/dynsyslab-public` installs the
  dependency-free `dynsyslab` client that drives the compiled CLI.
- **First calculation:**
  `dynsyslab-cli doctor` then
  `dynsyslab-cli run <example> --archive-root <a new folder>`, where `<example>` is
  `resources\examples\duffing-free-solve.json` in an installed release (under
  `%LOCALAPPDATA%\Programs\DynSysLab`) or `examples/native/duffing-free-solve.json` in this checkout.
- **Build from source:** [BUILD_NATIVE.md](BUILD_NATIVE.md).

## Where to read next

| Goal | Read |
|---|---|
| The whole suite at a glance, with diagrams | [OVERVIEW.md](OVERVIEW.md) |
| A five-minute workflow and the agent checklist | [START_HERE.md](START_HERE.md) |
| Package boundaries and engineering rules | [ARCHITECTURE.md](../ARCHITECTURE.md), [AGENTS.md](../AGENTS.md) |
| Jobs and pipelines | [architecture/PIPELINE_JOBS.md](architecture/PIPELINE_JOBS.md) |
| Estimation, control, linear systems | [architecture/ESTIMATOR_MODULE.md](architecture/ESTIMATOR_MODULE.md), [architecture/LINEAR_SYSTEMS_TOOLKIT.md](architecture/LINEAR_SYSTEMS_TOOLKIT.md) |
| Figures | [architecture/VIEWER_CLI.md](architecture/VIEWER_CLI.md) |

## Things you can ask your agent

- "Install DynSysLab and run the Duffing example."
- "Simulate the Lorenz system and show me its phase portrait."
- "Fit a Koopman model to these trajectories and check how well it predicts."
- "Design an LQR controller for this K-ROM and show me its Nyquist plot."
