# RPE echo-only KLQR: ψ(e) against ψ(ξ₁, ξ₂) on three targets

This pack regenerates, from nothing, a comparison of two ways to control the radius of a viscous
Rayleigh–Plesset bubble from its far-field acoustic echo alone. Both use KLQR on a Koopman
reduced-order model (K-ROM); they differ only in the eigenfunction the controller uses.

- **ψ(e):** eigenfunctions of the echo alone, evaluated directly from the measured echo, with no
  estimator.
- **ψ(ξ₁, ξ₂):** eigenfunctions of the first two causal echo-history coordinates. ξ is a frozen
  Hankel-SVD projection of the last 50 echo samples (decimation 10). The eigenvalues come from
  H-DMD (Arbabi & Mezić 2017) and the eigenfunctions from SINDy (Kaiser, Kutz & Brunton 2021).

Each is run on three radius targets:

| Target | Radius r(τ) |
|---|---|
| Mild | 1 + 0.05 sin(1.893 τ) |
| Strong | 1.5 + 0.4 sin(1.893 τ), ramped in from r = 1 over 10 τ |
| Hard | 1.12 plus four tones (0.95, 2.46, 3.97 and 5.49 rad per τ) |

## Results it reproduces

Radius error, the RMSE as a percentage of the target's own oscillation, over τ = 0 to 300:

| Target | ψ(e) | ψ(ξ₁, ξ₂) |
|---|---|---|
| Mild | 36.6 %, bounded; troughs flatten where ψ(e) folds (e ≈ 0.155) | 8.4 % |
| Strong | diverges at τ ≈ 3.1 | 17.6 % |
| Hard | diverges at τ ≈ 0.005 | 27.6 % |

A function of the echo alone cannot separate the two passes of each oscillation through the same
echo value, so ψ(e) folds; two history coordinates trace the circle an oscillatory eigenfunction
needs.

## Figures

For each run, Viewer renders:

- `fig-r-*`: the radius and its target;
- `spec-r-*`: their amplitude spectra;
- `fig-e-*`: the echo and its target, both lifted offline through (r, r′) by `control_target`;
- `spec-e-*`: their spectra;
- `fig-psi-*`: the controller's eigenfunction (mode 0) along the run and along its lifted target.
  Re ψ is plotted on the left axis and Im ψ on a stacked right axis.

ψ(e) on the hard target diverges too quickly to have a spectrum.

## Run it

See [docs/DATA_PACKS.md](../../docs/DATA_PACKS.md). It takes about 12 minutes on one machine.

```powershell
& $cli datapack validate .
& $cli run jobs\00-reproduce-all.job.json --archive-root D:\pack-run
```

## Layout

- `jobs/10-psi-e/`: the ψ(e) experiment.
- `jobs/20-psi-xi/`: the ψ(ξ₁, ξ₂) experiment.
- `jobs/00-reproduce-all.job.json`: both experiments, plus every figure stage: offline lifts,
  `figure_dataset` tables and `viewer_figure` renders.
