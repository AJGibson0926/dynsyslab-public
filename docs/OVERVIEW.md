# DynSysLab overview

DynSysLab is a native C++ suite for simulating, analysing, modelling and controlling dynamical
systems. One versioned archive holds every trajectory and every derived object, and three
applications work on it: the **Solver / Post-processor** (`DynSysLab.exe`) runs jobs, the
**Viewer** (`DynSysLabViewer.exe`) displays results, and the command-line tools
(`dynsyslab-cli`, `DynSysLabViewerCli`) do both headlessly for scripts and AI agents.

A typical study goes from equations to a controller in one reproducible chain:

1. **Simulate** a system from a library of ODEs, delay, lattice and PDE models, free or driven.
2. **Analyse** the trajectories: decompositions, Fourier spectra, phase space, Lyapunov screens.
3. **Identify a model**: DMD-family, Hankel/H-DMD, SINDy and spectral Koopman methods produce a
   Koopman reduced-order model (KROM) — eigenvalues, eigenfunctions, readouts and input maps.
4. **Design estimation and control** on that model: Koopman LQE (KLQE) estimates the lifted
   state from sensors, Koopman LQR/MPC (KLQR/KMPC) drives it to a target; classical LQR, PID,
   integral action and a linear-systems toolkit (gramians, balanced truncation, LQG) sit beside.
5. **Close the loop** in a new solve, archive everything, and compare runs in the Viewer.

Every step is a job (`dynsyslab.compute-job/v1`, JSON) that can be planned, validated, run and
re-run; derived archives point back to their sources, so each figure traces to its data.

---

## 1. Structure at a glance

```mermaid
flowchart LR
    subgraph APPS["Applications and tools"]
        SOLVERAPP["Solver / Post-processor<br/>DynSysLab.exe"]
        CLI["dynsyslab-cli<br/>validate · plan · run · queue"]
        PY["optional Python client"]
    end
    subgraph WF["Workflow service"]
        JOBS["Jobs<br/>solve · postprocess · pipeline"]
        PIPE["Experiments<br/>staged studies"]
        QUEUE["Durable queue"]
        HIER["Hierarchy<br/>Batch → System → Cluster → Trajectory"]
    end
    subgraph SCI["Scientific core (C++ libraries)"]
        SYS["Systems<br/>models · observables · targets"]
        SOLVE["Solvers<br/>RK45 · BDF · sampled control"]
        ANA["Analysis & identification<br/>decompositions · DMD · SINDy · Koopman"]
        KROM["Koopman models (KROM)"]
        EST["Estimation<br/>KLQE · observability"]
        CTRL["Control<br/>KLQR · KMPC · LQR · PID · LQI"]
        LIN["Linear-systems toolkit<br/>gramians · balancing · LQG"]
        NUM["Numerics<br/>linear algebra · Riccati · FFT"]
    end
    subgraph OUT["Data and display"]
        ARCH[("Archives<br/>versioned NPZ + index")]
        VIEWERAPP["Viewer<br/>DynSysLabViewer.exe"]
        VCLI["DynSysLabViewerCli<br/>headless figures"]
    end
    EXT["Extension DLLs<br/>systems · controllers · methods"]

    SOLVERAPP --> JOBS
    CLI --> JOBS
    PY --> CLI
    JOBS --> PIPE
    JOBS --> QUEUE
    JOBS --> HIER
    JOBS --> SOLVE
    JOBS --> ANA
    SYS --> SOLVE
    ANA --> KROM
    KROM --> EST
    KROM --> CTRL
    LIN --> EST
    LIN --> CTRL
    EST --> SOLVE
    CTRL --> SOLVE
    NUM --> LIN
    SOLVE --> ARCH
    ANA --> ARCH
    ARCH --> ANA
    ARCH --> VIEWERAPP
    ARCH --> VCLI
    EXT -. C ABI .-> SYS
    EXT -. C ABI .-> CTRL
```

The dependency direction is one-way: applications compose the workflow service, the workflow
composes scientific libraries, and only producers write archives. The Viewer is read-only.
Systems, controllers and post-processing methods can be added as separately compiled DLLs
through a versioned C ABI, without rebuilding the core.

---

## 2. Systems and solvers

```mermaid
flowchart LR
    subgraph LIB["System library"]
        ODE["ODE systems<br/>RPE · Marmottant · Duffing · Lorenz<br/>pendula · Hodgkin–Huxley · Kuramoto ..."]
        OTHER["Delay · lattice (Ising/Potts) · PDE"]
        OBS["Observable catalogue<br/>e.g. far-field acoustic echo,<br/>energies, pressures"]
        DRV["Drivers & targets<br/>waveforms · pulses · ramps"]
    end
    IC["Initial conditions<br/>explicit · grids · clusters"]
    CTL["Control policy<br/>continuous or sampled clock"]
    INT["Integrators<br/>Euler · RK4 · Dormand–Prince · BDF"]
    OUT[("Solution archive<br/>time · states · controls · targets<br/>· estimator traces · metadata")]

    ODE --> INT
    OTHER --> INT
    IC --> INT
    DRV --> CTL
    OBS --> CTL
    CTL <--> INT
    INT --> OUT
```

- A system supplies its right-hand side, Jacobian, input map, validity domain and an
  observable catalogue with analytic derivatives. Coordinates (dimensional or dimensionless)
  are explicit and recorded.
- Solvers integrate many trajectories per job. A controller is sampled on its own clock;
  with `solver.on_failure: "truncate"` a run that leaves the validity domain is archived up
  to the failure, with the reason.
- Campaign axes (parameter sweeps, initial-condition families) are part of the same job
  contract and run directly or through the durable queue.

---

## 3. Workflow: jobs, experiments and the hierarchy

```mermaid
flowchart TB
    JSON["Job JSON<br/>dynsyslab.compute-job/v1"] --> VAL["validate / plan<br/>no computation"]
    VAL --> RUN["run"]
    RUN --> K1["solve<br/>trajectories"]
    RUN --> K2["postprocess<br/>methods on sources"]
    RUN --> K3["pipeline<br/>stages with references"]
    K3 -->|"${stage.archive}"| K1
    K3 -->|"${stage.archive}"| K2
    K1 --> A[("Archives")]
    K2 --> A
    A --> M["Provenance manifest<br/>parameters · sources · versions"]
    subgraph H["Archive hierarchy"]
        B["Batch"] --> S["System + parameters"] --> C["Cluster"] --> T["Trajectory"] --> TS["TimeSeries"]
    end
    A --- H
```

- **Jobs** are deterministic and versioned. `plan` shows exactly what will run.
- **Experiments** chain stages (e.g. solve → H-DMD → KROM extraction → mode selection → closed
  loop); a stage refers to an earlier stage's archive, parameters can be overridden from the
  command line, and identical stages are reused.
- **Data roles** (free/forced, training, identification, reference) are metadata on archives,
  checked when a method requires, for example, unforced data for eigenfunctions.

---

## 4. Analysis and model identification

```mermaid
flowchart LR
    SRC[("Trajectories")] --> DEC["Decompositions<br/>SVD · PCA · POD · robust/bagged PCA"]
    SRC --> HANK["Hankel / delay embeddings<br/>Hankel SVD/POD · HAVOK"]
    SRC --> DMD["DMD family<br/>DMD · DMDc · EDMD(c) · H-DMD<br/>ResDMD · mpEDMD · TLS/Hermitian/mult DMD"]
    SRC --> SINDY["SINDy<br/>dynamics · eigenfunctions<br/>coordinate transforms"]
    SRC --> SPEC["Spectral Koopman<br/>measures · RKHS · stochastic"]
    DMD --> KX["KROM extraction<br/>mode selection by energy/rules"]
    HANK --> KX
    SINDY --> KX
    SPEC --> KX
    KX --> VALID["Validation<br/>semigroup / invariance tests"]
    KX --> MS["mode_selection<br/>eigenfunction subsets for<br/>estimation + control"]
    KX --> IM["Input maps B_ψ<br/>regressed from forced data"]
    VALID --> KROM[("KROM archive")]
    MS --> KROM
    IM --> KROM
```

- Methods run as post-processing stages on archived trajectories and write derived archives
  that reference their sources.
- A **KROM** stores continuous eigenvalues, eigenfunctions (polynomial, analytic-coordinate or
  scattered/interpolated values), named linear **readouts** (e.g. state, measurement) and
  **control maps** (constant, bilinear or exact B_ψ).
- `mode_selection` scores subsets of eigenfunctions with gramians, Hankel singular values and
  invariance tests, so an estimator and a controller can be built from the same model.

---

## 5. The Koopman model (KROM)

```mermaid
flowchart LR
    X["state x"] -->|"eigenfunctions ψ(x)<br/>(or of measured coordinates)"| Z["lifted state z"]
    Z -->|"ż = Λ z + B_ψ u"| Z
    Z -->|"readout C_y"| Y["measurement y ≈ y₀ + C_y z"]
    Z -->|"readout C_x"| XH["state x ≈ x₀ + C_x z"]
    U["input u"] -->|"B_ψ: constant · bilinear · exact"| Z
```

The KROM linearizes the dynamics in eigenfunction coordinates. The same object drives
prediction, estimation (through C_y) and control (through Λ and B_ψ).

---

## 6. Estimation and control

```mermaid
flowchart LR
    PLANT["Plant<br/>(simulated system)"] -->|"sensors only"| KLQE["KLQE<br/>ẑ̇ = Λẑ + B_ψu + L(y − ŷ)"]
    KLQE -->|"ẑ (lifted estimate)"| KLQR["KLQR / KMPC<br/>u = −K(ẑ − ψ_ref)"]
    KLQE -->|"x̂ = x₀ + C_x ẑ"| SF["any state-feedback controller<br/>(feedback: estimated_state)"]
    KLQR --> AUG["optional: integral action (LQI)<br/>· PID correction · cap · ramp"]
    SF --> AUG
    AUG -->|"u"| PLANT
    TGT["Target<br/>state · observable · lifted"] --> KLQR
    TGT --> SF
    LIN["Linear-systems toolkit<br/>gramians · HSV · balanced truncation<br/>pole placement · LQR/LQE/LQG · H2/H∞"] -.-> KLQE
    LIN -.-> KLQR
```

- **Controllers**: prescribed waveforms and pulses, PID, LQR/SDRE, MPC, KLQR and KMPC; all
  share one solver interface and an optional actuator cap.
- **Estimation**: KLQE uses only the configured sensors; it never evaluates eigenfunctions at
  the true state. Output feedback can share the controller's KROM (`feedback: estimate`) or
  hand a decoded state x̂ to any state-feedback controller (`feedback: estimated_state`).
- **Integral action and PID** are general components with KLQR adapters on any KROM readout.
- Analysis methods (`observability`, `controllability`, `linear_systems`) score models
  before they are used in a loop.

---

## 7. Archives and the Viewer

```mermaid
flowchart LR
    PROD["Producers<br/>solves · post-processing · experiments"] -->|"write once"| AR[("Archives<br/>NPZ arrays + metadata<br/>SQLite index")]
    AR -->|"read-only"| V["Viewer<br/>data tree · panels · comparisons"]
    AR -->|"read-only"| VC["Viewer CLI<br/>render · inspect · time-series"]
    V --> FIG["Figures<br/>time series · phase space · spectra<br/>control · estimator · KROM views"]
    VC --> FIG
    V -.-> CACHE["Disposable display cache"]
```

- Completed archives are immutable; derived results reference them instead of copying data.
- The Viewer reads the same archive tree; its CLI renders any panel headlessly, with time
  windows, dimensional units, log scales, axis ranges and multi-run comparisons.

---

## 8. Example: analyse and control a system

```mermaid
flowchart LR
    S1["1 · Solve free decays<br/>(training data)"] --> S2["2 · H-DMD / SINDy<br/>eigenvalues + eigenfunctions"]
    S2 --> S3["3 · KROM extraction<br/>+ mode selection"]
    S4["4 · Forced runs<br/>(identification data)"] --> S5["5 · Regress B_ψ<br/>and readouts"]
    S3 --> S5
    S5 --> S6["6 · Closed-loop solve<br/>KLQE → KLQR on a target"]
    S6 --> S7["7 · Compare in the Viewer<br/>vs exact inverse / full-state"]
```

```powershell
dynsyslab-cli run examples\native\experiments\rpe-echo-mode-selection\range-pool.json --archive-root archives
DynSysLabViewerCli render --source archives\solutions\rpe\<run>.npz --view "r versus time" --output radius.png
```

Worked case studies (bubble dynamics with echo-only sensing, exact linear references,
coordinate and history models) live under `examples/native/case-studies/` and
`docs/reports/`, each with the jobs that reproduce its figures.

---

## Where to go next

| Goal | Read |
|---|---|
| Install and first solve | [INSTALL.md](INSTALL.md), [START_HERE.md](START_HERE.md) |
| Write jobs and experiments | [architecture/EXPERIMENTS.md](architecture/EXPERIMENTS.md) |
| Koopman models | [architecture/EIGENFUNCTION_MODEL_API.md](architecture/EIGENFUNCTION_MODEL_API.md), [architecture/MODE_SELECTION.md](architecture/MODE_SELECTION.md) |
| Estimation and control | [architecture/ESTIMATOR_MODULE.md](architecture/ESTIMATOR_MODULE.md), [architecture/INTEGRAL_AND_PID_CONTROL.md](architecture/INTEGRAL_AND_PID_CONTROL.md), [architecture/LINEAR_SYSTEMS_TOOLKIT.md](architecture/LINEAR_SYSTEMS_TOOLKIT.md) |
| Figures | [architecture/VIEWER_CLI.md](architecture/VIEWER_CLI.md) |
| Extending the suite | [BUILD_NATIVE.md](BUILD_NATIVE.md) |
| Package boundaries | [../ARCHITECTURE.md](../ARCHITECTURE.md) |
