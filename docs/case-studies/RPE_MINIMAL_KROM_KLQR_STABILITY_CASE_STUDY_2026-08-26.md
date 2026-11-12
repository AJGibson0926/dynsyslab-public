# RPE minimal-KROM KLQR sampled-control stability case study

**Case date:** 2026-08-26  
**Status:** preserved baseline; diagnosis complete; baseline artifacts must not be overwritten  
**System:** controlled Rayleigh--Plesset equation (RPE)  
**Controller:** state-dependent Koopman LQR (KLQR/SDRE-style)  
**ROM:** one nonredundant complex Koopman eigenfunction plus its conjugate dual  
**Principal finding:** a fast continuous-time KLQR pole becomes a negative, lightly damped sampled pole under the nominal 0.005-`tau` zero-order hold. The visible control “squiggle” near `tau = 15.9` is sampled high-gain chatter, not pressure saturation and not a Viewer artifact. A separate floating-point clock defect can cause duplicate or half-step-late Riccati updates; it is a secondary implementation artifact that must also be preserved when reproducing the original run.

This document is the authoritative human-readable record of the case. The companion [reproduction notes](../notes/RPE_MINIMAL_KROM_KLQR_STABILITY_REPRO_NOTES_2026-08-26.md) provide a shorter laboratory checklist.

## 1. Scientific question and conclusion

The controlled solution tracks a small two-tone RPE target well, but the pressure command develops a conspicuous alternating component just after `tau = 15`. The investigation asked whether that behavior came from:

1. the nonlinear RPE plant;
2. the extracted KROM or its state-dependent input map;
3. Riccati tuning;
4. pressure saturation;
5. integration or controller sampling; or
6. visualization.

The dominant mechanism is the interaction between aggressive continuous-time Riccati tuning and the explicit controller sample-and-hold clock:

- the controller uses `Q = I_2` and `R = [[1e-16]]`;
- the fastest frozen continuous-time closed-loop pole at `tau ~= 15.93` is `-342.5695` per dimensionless time unit;
- its continuous-time time constant is approximately `0.002919 tau`;
- the intended controller interval is `0.005 tau`, with the command held between updates;
- exact zero-order-hold discretization of the frozen local model gives the corresponding sampled pole `-0.7362623`;
- the negative sign reverses the fast error component at every controller update, producing the observed alternating pressure command.

The sampled local model remains stable at the representative point because the pole magnitude is below one. Therefore “unstable-looking” or “chattering” is more precise than “locally unstable” for the preserved baseline. A tenfold increase in `Q` at unchanged `R` does make the corresponding sampled pole genuinely unstable (`-4.49093` at the same point).

## 2. Immutable artifact inventory

Paths below are relative to the DynSysLab source root. SHA-256 hashes are the primary identity checks.

| Role | Immutable artifact | Bytes | SHA-256 |
|---|---|---:|---|
| Controlled baseline | `archives/solutions/rayleigh-plesset/rpe-equilibrium-minimal-krom-klqr-multisine-tau25__20260826T170659.584050Z.npz` | 1,043,506 | `91a1b9255dcde946846712cb530ffd5ef80bc88877d410c5ed32e302aef43fdb` |
| Extracted minimal KROM | `archives/derived/koopman-krom-extract/rpe-31-free-natural-period-dtau0p001-tau25-c000-t030__20260821T203315.886026Z__koopman-krom-extract__20260826T164339.807992Z__ebb5aedc.npz` | 1,173,124 | `16dbd0675707021aa94c6452ff0d978a04a4f6c423e0a204993bbee628363712` |
| Full SINDy KROM | `archives/derived/koopman-with-sindy/rpe-31-free-natural-period-dtau0p001-tau25-c000-t030__20260821T203315.886026Z__koopman-with-sindy__20260821T222933.502003Z__50b83a95.npz` | 419,158,070 | `ac4d3e6e4062a417aee538fdcbf089f84d0238a0e3c850e79d75bef8f1e55e2b` |
| Full H-DMD result | `archives/derived/koopman-hankel-dmd/rpe-31-free-natural-period-dtau0p001-tau25-c000-t030__20260821T203315.886026Z__koopman-hankel-dmd__20260821T222910.100953Z__d155f35a.npz` | 378,497,096 | `bfa7c74db990cfe19abe420fa1370e9702133204fa189bfe5bfd6f7ffc61632a` |
| Free training trajectory | `archives/solutions/rayleigh-plesset/rpe-31-free-natural-period-dtau0p001-tau25-c000-t030__20260821T203315.886026Z.npz` | 342,638 | `92432b83328e9d83a6b7822b541c2d844106286ca76a34c901cc8bf089450090` |

The controlled archive contains:

| Array | Shape | Meaning |
|---|---:|---|
| `values` | `(25001, 2)` | canonical dimensional state samples `(R, dR/dt)` |
| `time` | `(25001,)` | dimensional seconds corresponding to `tau = 0:0.001:25` |
| `control_values` | `(25001, 1)` | held pressure command in Pa, aligned with `time` |
| `control_call_time` | `(100000,)` | all RK4 RHS-call times, including repeated stage times |
| `control_call_values` | `(100000, 1)` | pressure returned to every RHS call |
| `control_target_values` | `(25001, 2)` | target in canonical dimensional state coordinates |

The archive is the ground truth. Do not replace it when running tuning experiments; use distinct labels and immutable derivative archives.

## 3. Code and runtime identity

The repository was on branch `main` at commit:

```text
7037594f4db6874d3469a032b0cf245e446dac81
```

The working tree was dirty. Consequently, the commit alone is not sufficient to recreate the case. Relevant file hashes at documentation time were:

| Source module | SHA-256 |
|---|---|
| `dynsyslib/lib/rpe.py` | `5c298bae8bb0a3287e4c06be1c9b039b9800bb1cf30cdbbec377a7269e8fc259` |
| `control/koopman.py` | `7845a30847c6c8d459ebfbe443d4d8984425781855cb9fb2e52f652978012a73` |
| `control/targets.py` | `a60c20e3509447571d314d56cc24ed821bd1944557f411dbe1be7098d4da1be2` |
| `solver/ode.py` | `5248c364de3750681290184ec5a9dc5ec7c0a47b39a7015d4013dd5695790d96` |
| `rom/koopman.py` | `dec914821f77ff7397e00d7d551f1e1aadfd1e39b8bf1dd08c770800f47ceb6d` |
| `rom/storage.py` | `6cf6c9ac3b6b34f96d941a271ffbf86a1911a7f91cd25c8eb629f598521ef451` |
| `postprocessing/koopman/pipeline.py` | `13e635cfcc454712c6346a484d22bece86de4119c1edb1b4f98bd03c581ae84a` |
| `postprocessing/koopman/koopman_sindy.py` | `d2dca4459d10dedd60d45d2086dfa30d22674aeacfc5251706645cf6820276f9` |
| `postprocessing/koopman/dmd.py` | `62a8e59fab525d2fd0c32fee94a2084788fb46571f9ff798864ec7b158303428` |

Runtime:

```text
Python 3.12.1 (Windows, MSC v.1937 64 bit)
NumPy 2.5.1
SciPy 1.18.0
Solver execution backend: python-reference
```

## 4. RPE dynamical system

### 4.1 Physical parameters

The case uses the `air_water_20c` preset.

| Symbol | Code field | Value | Unit |
|---|---|---:|---|
| `R0` | `equilibrium_radius_m` | `1.0e-5` | m |
| `kappa` | `preset.kappa` | `1.3333` | unitless |
| `p0` | `ambient_pressure_pa` | `100000` | Pa |
| `pV` | `preset.vapor_pressure_pa` | `2330` | Pa |
| `rhoL` | `preset.rho_kg_m3` | `998` | kg/m^3 |
| `sigma` | `preset.sigma_n_m` | `0.0725` | N/m |
| `muL` | `preset.mu_pa_s` | `0.001` | Pa s |
| `pG0` | derived equilibrium gas pressure | `112170` | Pa |

The equilibrium gas pressure is

```text
pG0 = p0 - pV + 2 sigma/R0 = 112170 Pa.
```

There are no `CL` or `CG` parameters in this incompressible RPE model.

### 4.2 Controlled dimensional equation

With radius `R`, radial velocity `V = dR/dt`, and controller pressure `p_c`, the implementation is

```text
dR/dt = V
dV/dt = [pG0 (R0/R)^(3 kappa) + pV - p0
         - 2 sigma/R - 4 muL V/R] / (rhoL R)
         - 3 V^2/(2 R)
         + p_c/(rhoL R).
```

The control-affine input matrix is

```text
B(R) = [0, 1/(rhoL R)]^T.
```

Positive `p_c` therefore enters with the positive sign shown above. The uncontrolled driver is zero for this case.

## 5. Nondimensionalization

The convention is `rpe-natural-period/v1`:

```text
r   = R/R0
tau = t/T0 = f0 t
u   = dr/dtau = (T0/R0) dR/dt.
```

The requested natural frequency and characteristic time are

```text
f0 = [1/(2 pi R0 sqrt(rhoL))] sqrt(3 kappa pG0 - 2 sigma/R0)
T0 = 1/f0.
```

For this case:

| Quantity | Value |
|---|---:|
| `f0` | `331958.7786861145 Hz` |
| `T0` | `3.0124222168727627e-6 s` |
| `P0 = rhoL (R0/T0)^2` | `10997.623748528325 Pa` |

The dimensionless groups stored in both training and controlled archives are:

| Group | Definition | Value |
|---|---|---:|
| gas pressure | `pG0/P0` | `10.199476047269785` |
| vapor pressure | `pV/P0` | `0.21186394927465987` |
| ambient pressure | `p0/P0` | `9.092873359427463` |
| capillarity | `2 sigma/(R0 P0)` | `1.318466637116982` |
| viscosity | `4 muL/(T0 P0)` | `0.12073836540572193` |

Define `p_c* = p_c/P0`. The integrated dimensionless system is

```text
dr/dtau = u
du/dtau = [Pi_g r^(-3 kappa) + Pi_v - Pi_0
           - Pi_sigma/r - Pi_mu u/r - 1.5 u^2]/r
           + p_c*/r.
```

The solver integrates this form but publishes dimensional `time` and `values` in the archive. Analysis and Viewer recover `tau`, `r`, `u`, and `p_c*` from the stored nondimensionalization metadata.

## 6. Free trajectory used to derive the KROM

The minimal KROM ultimately comes from one free-response trajectory, not from the controlled trajectory:

```text
label: rpe-31-free-natural-period-dtau0p001-tau25-c000-t030
initial physical state: (2.3334e-5 m, 0 m/s)
initial dimensionless state: (2.3334, 0)
interval: tau in [0, 25]
published spacing: delta_tau = 0.001
published samples: 25001
uncontrolled/free response: yes
```

The source solve used dimensionless RPE equations with `DOP853`, `rtol = 1e-9`, and `atol = 1e-12`. It produced 1,046 raw adaptive samples and was resampled to the uniform 25,001-sample analysis series.

## 7. H-DMD derivation

The full H-DMD stage used the standard trajectory policy documented in [HDMD_SINDY_PIPELINE.md](../HDMD_SINDY_PIPELINE.md):

| Setting | Value |
|---|---:|
| observable | dimensionless radius, state index 0 |
| native/analysis spacing | `delta_tau = 0.001` |
| analysis samples | `22000` |
| sampling window | `20000` samples (`19.999 tau`) |
| matrix size / delays | `2000` |
| Hankel `X`, `Y` shape | `20000 x 2000` each |
| window duration per delay vector | `1.999 tau` |
| delay stride | 1 sample |
| shifted pair | one-sample shift |
| SVD driver | LAPACK `gesvd` |
| rank policy | absolute singular-value tolerance |
| tolerance | `1e-10` |
| retained rank | `2000` (full spectrum; no 12-mode truncation) |
| exact modes | false; projected modes `U W` |

The H-DMD operator produces unitless discrete eigenvalues `mu`. Continuous generator eigenvalues are obtained only through

```text
lambda = log(mu)/delta_tau.
```

The recorded round-trip error `|exp(lambda delta_tau) - mu|` over the full transferred spectrum was at most `7.108895957933346e-16`.

Recorded H-DMD timing:

| Stage | Seconds |
|---|---:|
| SVD | `87.1963751` |
| projected operator | `4.2905577` |
| eigendecomposition | `36.5186477` |
| mode construction | `13.1435694` |
| total | `141.1491682` |

The H-DMD archive contains 2,000 eigenvalues, a `2000 x 2000` real projected operator, `20000 x 2000` complex modes, and 2,000 initial modal amplitudes.

## 8. SINDy eigenfunction derivation

SINDy did not create a new spectrum. It received the exact 2,000 H-DMD discrete eigenvalues, converted each to a continuous generator eigenvalue with the logarithmic map above, and fit polynomial Koopman eigenfunctions for those supplied generators.

| Setting | Value |
|---|---:|
| method | `koopman-with-sindy` |
| input trajectories | 1 |
| state coordinates | dimensionless `(r, u)` |
| polynomial degree | 4 |
| feature count | 15 |
| constant included | yes |
| state normalization | yes |
| sparse threshold | `0.05` |
| maximum iterations | `20` |
| null-space backend | guarded Gram-Hermitian eigensolve (`gram-eigh`) |
| generator data | dimensionless system RHS |
| conjugate reuse | 999 fits reused; 1001 direct fits |
| full eigenfunctions | 2,000 |
| nonredundant KROM coordinates | 1,001 |

The normalized variables used by every archived polynomial are

```text
z1 = (r - 1.4549767945605923) / 0.43850754420434557
z2 = (u + 0.0372157023019093) / 3.348671331285554.
```

The feature order is

```text
1, z1, z2,
z1^2, z1 z2, z2^2,
z1^3, z1^2 z2, z1 z2^2, z2^3,
z1^4, z1^3 z2, z1^2 z2^2, z1 z2^3, z2^4.
```

Recorded SINDy construction time was `1.4941908 s`, including `0.4343427 s` for sparse null-space fits and `0.1154150 s` for batched eigenfunction evaluation. This timing describes the already-created archive, not Viewer display-time fitting.

## 9. Minimal KROM extraction

The reserved archive label is:

```text
rpe-natural-fundamental-single-eigenfunction
```

Extraction requested full-spectrum source index `433` and included its conjugate dual at index `434`. The archive therefore stores two eigenvalues and two exactly conjugate coefficient rows, while the KROM exposes one nonredundant complex eigenfunction as two real coordinates.

### 9.1 Eigenvalue pair

| Domain | Positive-frequency value | Conjugate |
|---|---|---|
| discrete | `0.9998540377155265 + 0.005598164269206836 i` | complex conjugate |
| continuous, per `tau` | `-0.13029888671192555 + 5.598923003706494 i` | complex conjugate |

Derived descriptors:

| Descriptor | Value |
|---|---:|
| `|mu|` | `0.9998697096018193` |
| angle of `mu` | `0.005598923003706493 rad/sample` |
| oscillation frequency | `0.8910962720308107 cycles/tau` |
| oscillation period | `1.1222132011853192 tau` |
| decay time | `7.67466265625795 tau` |
| modal Q-factor | `21.484922645904906` |

### 9.2 Exact polynomial eigenfunction

For the positive-frequency eigenvalue, `psi(r,u) = Theta(z1,z2) c` with:

| Feature | Real coefficient | Imaginary coefficient |
|---|---:|---:|
| `1` | 0 | 0 |
| `z1` | 0 | 0 |
| `z2` | 0 | 0 |
| `z1^2` | 0 | 0 |
| `z1 z2` | `-0.28107901541197644` | 0 |
| `z2^2` | `0.03856723735576586` | `0.7299779079221112` |
| `z1^3` | 0 | 0 |
| `z1^2 z2` | `-0.20514884094488225` | `-0.0007860010228844603` |
| `z1 z2^2` | `0.030314223833558654` | `0.5739463259918604` |
| `z2^3` | 0 | 0 |
| `z1^4` | 0 | 0 |
| `z1^3 z2` | `-0.03772727235859941` | `-0.00027020988777596146` |
| `z1^2 z2^2` | `0.005963971428232476` | `0.11296202092036441` |
| `z1 z2^3` | 0 | 0 |
| `z2^4` | 0 | 0 |

The negative-frequency coefficient vector is the exact complex conjugate; the maximum archived conjugacy discrepancy is zero.

### 9.3 KROM representation and limitations

The nonredundant split coordinate is

```text
xi = [Re psi, Im psi]^T.
```

Its autonomous dimensionless generator is

```text
A_tau = [[-0.13029888671192555, -5.598923003706494],
         [ 5.598923003706494, -0.13029888671192555]].
```

The KROM is eigenfunction-only:

- lifted evolution: available;
- analytic control coupling: available;
- state reconstruction: unavailable;
- control representation: analytic state-dependent `B_psi(x)`;
- no learned linear or bilinear `B1/B2` matrices are used.

Analytic control coupling is

```text
B_psi(x) = grad(psi)(x) B(x),
```

split into real and imaginary rows. Because the polynomial gradient and physical RPE input matrix both depend on state, `B_psi` and the Riccati gain are recomputed on the controller clock.

## 10. Control target

The canonical target is `rpe_multisine_equilibrium_small`:

```text
r_d(tau) = 1 + 0.06 sin(4 tau) + 0.08 sin(3 tau + pi)
         = 1 + 0.06 sin(4 tau) - 0.08 sin(3 tau)

u_d(tau) = 0.24 cos(4 tau) + 0.24 cos(3 tau + pi)
         = 0.24 cos(4 tau) - 0.24 cos(3 tau).
```

It begins exactly at equilibrium:

```text
[r_d(0), u_d(0)] = [1, 0].
```

The controlled RPE solve is independently gated to require the physical initial state `(R0, 0)`. Thus the controller must push the bubble away from equilibrium; the solve cannot begin on a later target state.

For KLQR, the physical state target is first evaluated in dimensionless model time, then lifted through the same polynomial eigenfunction:

```text
xi_d(tau) = [Re psi(r_d,u_d), Im psi(r_d,u_d)]^T.
```

## 11. KLQR law and tuning

At each controller update, the implementation forms:

```text
A = A_tau/T0
B = B_psi(x) expressed per physical second and per Pa
```

and solves the continuous algebraic Riccati equation

```text
A^T P + P A - P B R^-1 B^T P + Q = 0.
```

The gain and pressure command are

```text
K = R^-1 B^T P
p_c = -K (xi - xi_d).
```

Exact preserved tuning:

```text
Q = [[1, 0],
     [0, 1]]

R = [[1e-16]]

pressure limit = +/-200000 Pa.
```

`Q` penalizes the real and imaginary eigenfunction-coordinate error. It does **not** directly penalize radius and radial-velocity errors. `R` penalizes dimensional pressure squared. In dimensionless pressure coordinates, the equivalent scalar weight is

```text
R* = R P0^2 = 1.209477281141942e-8.
```

This mixed scaling is central to the very large gain. A future normalized controller should express the entire cost in declared dimensionless coordinates.

## 12. Controlled solve configuration

| Setting | Value |
|---|---:|
| initial physical state | `(1e-5 m, 0 m/s)` |
| initial dimensionless state | `(1, 0)` |
| interval | `tau in [0,25]` = `t in [0,7.531055542181906e-5 s]` |
| integration method | fixed-step RK4 |
| integration step | `delta_tau = 0.001` = `3.012422216872763e-9 s` |
| nominal controller hold | `delta_tau_c = 0.005` = `1.5062111084363815e-8 s` |
| control updates per target interval | 5 solution steps per held command |
| state samples | 25,001 |
| RK4 RHS/controller records | 100,000 |
| state estimator | identity/full-state feedback |
| execution backend | Python reference |
| controller active window | none; active for full solve |

The controller clock is explicit and independent of RHS evaluation density. It is intended to recompute when the integration time reaches the next 0.005-`tau` boundary and otherwise return the held command. RK4 still calls the controlled RHS four times per state step, which explains the 100,000 controller-call records; those are not 100,000 distinct Riccati updates.

The preserved implementation computes the next boundary with `floor(t/step)`. Floating-point accumulated time and the clock's comparison tolerance can therefore re-enter the same nominal tick or miss it until the next RK4 half-stage. This produced a small but real departure from the ideal schedule, documented in Section 15.1.

## 13. Baseline performance

All errors below are measured after converting the archived dimensional state and target back to `(r,u)`.

| Metric | Value |
|---|---:|
| radius RMSE | `0.006360083992258022` |
| radial-velocity RMSE | `0.007339480878497758` |
| maximum absolute radius error | `0.011419017629996864` |
| maximum absolute velocity error | `0.022717496907171505` |
| mean pressure | `-2721.0146152919783 Pa` |
| pressure RMS | `23775.90402176294 Pa` |
| peak absolute pressure | `68040.66214461008 Pa` |
| dimensionless pressure RMS | `2.161912842757925` |
| peak absolute dimensionless pressure | `6.1868512417252965` |
| saturated samples | `0` |

The ±200 kPa limiter is therefore not the source of the oscillation.

## 14. Representative local stability calculation

At the representative controller boundary near the largest alternating jumps:

```text
tau = 15.92999999999661
t   = 4.79878859147729e-5 s
```

the archived state, target, and error are:

| Quantity | Radius component | Velocity component |
|---|---:|---:|
| dimensionless state | `1.1043053887412562` | `0.33476114547142466` |
| dimensionless target | `1.0959888348400604` | `0.3400997180367939` |
| state minus target | `0.008316553901195789` | `-0.005338572565369237` |
| physical state | `1.1043053887412562e-5 m` | `1.1112690100225886 m/s` |
| physical target | `1.0959888348400605e-5 m` | `1.1289908703098603 m/s` |

The held pressure is `38154.56984837685 Pa`, or `p_c* = 3.469346717147202`.

The physical-time lifted generator is

```text
A = [[  -43253.8592970512, -1858611.6422679997],
     [ 1858611.6422679997,   -43253.8592970512]] 1/s.
```

The state-dependent lifted input vector and Riccati quantities are

```text
B_psi = [[0.9542494769444892],
         [0.6194563192815491]]

P = [[ 1.5685454150012215e-7, -2.3632607747880300e-7],
     [-2.3632607747880300e-7,  3.8598093895507450e-7]]

K = [[3.2846820775937274e7, 1.3584295935456130e8]].
```

Additional conditioning indicators:

```text
||B_psi||_2 = 1.1376810606433931
cond([B_psi, A B_psi]) = 1.8596182518800036e6.
```

The large controllability-matrix condition number indicates poor numerical scaling in the split eigenfunction representation, although the input channel is not singular at the event.

The frozen continuous-time closed-loop poles are:

| Coordinate | Pole 1 | Pole 2 |
|---|---:|---:|
| per second | `-1.8603986044410318e6` | `-1.137189502609702e8` |
| per `tau` | `-5.604306088257247` | `-342.5694922455953` |

Both continuous-time poles are stable. The fast time constant is only `1/342.5695 = 0.002919 tau`, shorter than the `0.005 tau` hold.

## 15. Exact sampled local model

For a frozen `(A,B,K)` and hold `h`, the diagnostic used the exact zero-order-hold matrices

```text
A_d = exp(A h)
B_d = integral_0^h exp(A s) B ds
```

and inspected

```text
A_d,cl = A_d - B_d K.
```

At `tau ~= 15.93` and `h = 0.005 T0`, its eigenvalues are

```text
z_cl = {0.9723761687799314, -0.7362623165043962}.
```

The slow pole is positive and close to one. The fast pole is stable in magnitude but negative, so its component alternates sign on each control update. The nonlinear target continually excites it, making the pressure oscillation persistent and visually prominent.

Across `14 <= tau <= 17`, using frozen local models at each 0.005-`tau` update:

```text
sampled spectral radius range: 0.9714322 to 0.9723794
most negative sampled pole observed: approximately -0.93788
fraction with spectral radius above one: 0
```

This is local, frozen-model evidence. It is not a global nonlinear stability proof.

### 15.1 Preserved controller-clock defect

The controller call log contains 100,000 RK4-stage records. Counting every change in the returned command gives 5,058 changes, rather than the nominal 5,000 interval transitions. The intervals between command changes are:

| Observed interval in `tau` | Count |
|---:|---:|
| `0.005` | 4,969 |
| `0.0045` | 29 |
| `0.0005` | 29 |
| `0.0` (same solver time) | 29 |
| `0.004999999992` | 1 |

The 29 zero-time transitions are genuine duplicate controller evaluations at the same floating-point time but different RK4 stage states. A typical sequence near `tau = 1.01` is:

```text
tau=1.0099999999999996: -17419.30 -> -17122.96 Pa
tau=1.0099999999999996: -17122.96 -> -17138.87 Pa
tau=1.0104999999999995: -17138.87 -> -17128.64 Pa
```

Near the main chatter event, the changes are instead evaluated approximately half a solver step after the nominal boundaries:

```text
15.9005, 15.9055, 15.9105, ... tau.
```

The effective interval there is still approximately `0.005 tau`, so the negative sampled-pole mechanism remains the primary explanation. The clock defect changes the exact update phase and occasionally inserts extra Riccati evaluations; it can amplify or perturb a high-gain case and prevents the original run from being interpreted as a perfectly periodic digital controller.

The correct future implementation should use an integer tick index, advance it monotonically, and compute the boundary as `tick*step`, never as `floor(t/step)+1` after a tolerant early trigger. A regression test should assert exactly one controller update per tick for every fixed-step RK stage pattern.

## 16. Observed alternating pressure sequence

The largest command jumps in the solution-aligned archive occur around `tau = 15.9`. The table reports adjacent stored solution samples; the underlying controller evaluations occur near the half-step times described above:

| Boundary | Before (Pa) | After (Pa) | Jump (Pa) | Dimensionless jump |
|---|---:|---:|---:|---:|
| `15.925 -> 15.926` | `8329.1162` | `38154.5698` | `+29825.4536` | `+2.71199` |
| `15.930 -> 15.931` | `38154.5698` | `8543.6140` | `-29610.9558` | `-2.69249` |
| `15.935 -> 15.936` | `8543.6140` | `38011.5749` | `+29467.9608` | `+2.67948` |
| `15.940 -> 15.941` | `38011.5749` | `9314.0759` | `-28697.4990` | `-2.60943` |
| `15.945 -> 15.946` | `9314.0759` | `37271.8881` | `+27957.8122` | `+2.54217` |

Each command is held for approximately five `delta_tau = 0.001` solution intervals in this region. The alternating pattern is synchronized to the controller clock, not to plotting pixels or saturation.

The target near `tau = 15.93` is

```text
r_d = 1.09598883484
u_d = 0.34009971803
d u_d/dtau = -1.18977922227.
```

This is a high-velocity, high-curvature portion of the two-tone target. The rapidly changing demand and small model/state error excite the high-gain sampled mode, explaining why the chatter becomes especially visible there rather than at the beginning of the run.

## 17. Tuning and sampling sensitivity

The following are frozen-model counterfactuals at `tau ~= 15.93`; they are not replacement nonlinear simulations.

| `Q` | `R` | Fast continuous pole per `tau` | Fast sampled pole at `delta_tau_c=0.005` | Interpretation |
|---|---:|---:|---:|---|
| `I` | `1e-16` | `-342.5695` | `-0.73626` | preserved baseline; stable alternating mode |
| `I` | `1e-15` | `-107.9291` | `+0.45238` | removes sign reversal; less aggressive |
| `I` | `1e-14` | `-32.7966` | `+0.83274` | smoother but slower control |
| `0.1 I` | `1e-16` | `-107.9291` | `+0.45238` | essentially same `Q/R` ratio as row above |
| `10 I` | `1e-16` | `-1083.7140` | `-4.49093` | sampled local instability |

The equivalence between scaling `Q` down by ten and scaling `R` up by ten follows from the relative cost scaling. Increasing `Q` or decreasing `R` makes control cheaper relative to state error, raises gain, and worsens sampled stability. Increasing `R` reduces actual actuator effort; decreasing it does the opposite.

At baseline `Q,R`, changing only the controller hold gives the following worst fast sampled pole over `14 <= tau <= 17`:

| Controller interval | Most negative sampled pole |
|---:|---:|
| `0.005 tau` | `-0.93788` |
| `0.002 tau` | `+0.23107` |
| `0.001 tau` | `+0.61658` |
| `0.0005 tau` | `+0.80855` |

Shortening the hold eliminates sign reversal while retaining the continuous-time cost. The apparent movement of the positive pole toward one as the step shrinks is expected for a fixed continuous decay rate measured over a shorter discrete interval.

Evaluating alternate gains along the preserved trajectory, rather than rerunning the nonlinear plant, gave the following local pressure-jump indicators over `15.75 <= tau <= 16.05`:

| `Q` | `R` | Maximum adjacent 0.005-`tau` command change on preserved path |
|---|---:|---:|
| `I` | `1e-16` | `17566 Pa` |
| `I` | `1e-15` | `5716 Pa` |
| `I` | `1e-14` | `1937 Pa` |
| `10 I` | `1e-16` | `55001 Pa` |

These counterfactual commands do not include the changed plant trajectory that a new gain would produce. They are sensitivity evidence only.

### 17.1 Completed nonlinear follow-up runs

Two full 25-`tau` nonlinear experiments were subsequently run with the same plant, equilibrium initial condition, target, minimal KROM, RK4 step, pressure bound, and Python backend. Both are new immutable archives; the baseline was not modified.

| Case | `R` | `delta_tau_c` | Archive SHA-256 | Solve time |
|---|---:|---:|---|---:|
| conservative combined | `1e-15` | `.001` | `3d31831d6835bc84af24975f447fbe1b22f04a5bf2124c74fcf1eada8c8b6339` | `104.75 s` |
| sampling-only / selected | `1e-16` | `.001` | `b94345fd39cba65f99c634937cc7093955b3724475982e0fd49bb64fb5af79ad` | `91.14 s` |

Artifact names:

```text
archives/solutions/rayleigh-plesset/
rpe-equilibrium-minimal-krom-klqr-multisine-tau25-R1e-15-dtauc0p001__20260826T175222.433924Z.npz

archives/solutions/rayleigh-plesset/
rpe-equilibrium-minimal-krom-klqr-multisine-tau25-R1e-16-dtauc0p001__20260826T175520.700735Z.npz
```

Measured comparison:

| Metric | Baseline `R=1e-16`, `.005` | Combined `R=1e-15`, `.001` | Sampling-only `R=1e-16`, `.001` |
|---|---:|---:|---:|
| radius RMSE | `0.00636008` | `0.01153073` | `0.00635676` |
| velocity RMSE | `0.00733948` | `0.02105565` | `0.00733403` |
| max radius error | `0.01141902` | `0.02353196` | `0.01140902` |
| max velocity error | `0.02271750` | `0.06381752` | `0.02268025` |
| pressure RMS (Pa) | `23775.90` | `22585.50` | `23762.90` |
| peak pressure (Pa) | `68040.66` | `61313.93` | `68792.98` |
| maximum adjacent command jump (Pa) | `29825.45` | `628.17` | `2919.34` |
| maximum jump, `14<=tau<=17` (Pa) | `29825.45` | `559.19` | `2544.84` |
| pressure total variation, `14<=tau<=17` (Pa) | `1046680.71` | `244508.48` | `264171.59` |
| saturated samples | 0 | 0 | 0 |

The conservative combined case removed nearly all chatter but degraded radius RMSE by 81.3% and velocity RMSE by 186.9%. It therefore did not satisfy the joint objective.

The sampling-only case satisfied both objectives:

- radius RMSE improved by `0.052%`;
- velocity RMSE improved by `0.074%`;
- maximum command jump fell by `90.21%`;
- maximum jump over `14<=tau<=17` fell by `91.47%`;
- control total variation over that interval fell by `74.76%`;
- pressure RMS fell by `0.055%` and peak pressure rose by `1.11%`.

At `tau ~= 15.93`, the sampling-only sampled poles were

```text
{0.99441183518778, 0.6573345727269536}.
```

The problematic pole is now positive, so the update-to-update sign reversal is removed. The slow pole is closer to one per sample only because each sample is five times shorter.

This experiment isolates the cause cleanly: increasing `R` reduced control authority and hurt tracking, while shortening the implemented control interval removed chatter without changing the continuous-time cost ratio. The selected follow-up configuration is therefore

```text
Q = I_2
R = [[1e-16]]
delta_tau_c = 0.001.
```

The floating clock defect remains. The two `.001` runs each recorded 25,296 command changes and 148 same-time changes, rather than an ideal 25,000 interval transitions. Thus the selected result improves substantially despite the defect; correcting the clock remains a separate required experiment.

Viewer artifact for the selected run:

```text
viewer_outputs/rpe-minimal-krom-klqr-multisine-controlled-dtauc0p001.html
```

### 17.2 Alternate single-eigenfunction experiment

A second eigenfunction was selected from the same original 2,000-mode H-DMD/SINDy KROM without refitting. The screening criteria were:

- negative continuous growth rate;
- positive-frequency representative with a conjugate dual;
- angular frequency closer to the target's `3` and `4 rad/tau` components;
- nonzero sparse polynomial coefficients;
- analytic `B_psi` of comparable magnitude to the index-433 mode.

Full-spectrum index `431` was the best match under those criteria. Index `432` is its conjugate. The original selected mode was index `433` plus conjugate `434`.

| Property | Original index 433 | Alternate index 431 |
|---|---:|---:|
| continuous `lambda` per `tau` | `-0.13029889 + 5.59892300 i` | `-0.21123031 + 4.82274786 i` |
| discrete `mu` | `0.99985404 + 0.00559816 i` | `0.99977717 + 0.00482171 i` |
| oscillation period (`tau`) | `1.12221` | `1.30282` |
| decay time (`tau`) | `7.67466` | `4.73417` |
| nonzero degree-4 polynomial coefficients | 6 | 6 |

The alternate pair was immutably reserved as:

```text
archives/derived/koopman-krom-extract/
rpe-31-free-natural-period-dtau0p001-tau25-c000-t030__20260821T203315.886026Z__koopman-krom-extract__20260826T180431.756269Z__f66c8d61.npz
SHA-256: f0cbfb21def5c4492b5c6a20b7f3a22ff224b41b3b80c57ed05645688d3d755c
reservation label: rpe-alternate-target-nearby-single-eigenfunction-index431
```

Everything else was held equal to the selected index-433 experiment: identical plant and parameters, target, equilibrium initial state, `Q=I`, `R=1e-16`, `delta_tau=delta_tau_c=.001`, RK4, pressure limit, backend, and 25-`tau` duration.

Alternate controlled archive:

```text
archives/solutions/rayleigh-plesset/
rpe-equilibrium-alt-krom-index431-klqr-multisine-tau25-R1e-16-dtauc0p001__20260826T180450.293836Z.npz
SHA-256: 925619324dbab3dab8143e28bcf93bf1de049c891b4d6363158e8f43e1a9125a
solve time: 119.81 s
```

Direct equal-setting comparison:

| Metric | Original index 433 | Alternate index 431 | Alternate change |
|---|---:|---:|---:|
| radius RMSE | `0.00635676` | `0.00660105` | `+3.84%` |
| velocity RMSE | `0.00733403` | `0.00790828` | `+7.83%` |
| max radius error | `0.01140902` | `0.01194708` | `+4.72%` |
| max velocity error | `0.02268025` | `0.02454989` | `+8.24%` |
| pressure RMS (Pa) | `23762.90` | `23689.72` | `-0.31%` |
| peak pressure (Pa) | `68792.98` | `68065.10` | `-1.06%` |
| maximum command jump (Pa) | `2919.34` | `2506.24` | `-14.15%` |
| max jump, `14<=tau<=17` (Pa) | `2544.84` | `1863.26` | `-26.78%` |
| pressure TV, `14<=tau<=17` (Pa) | `264171.59` | `256629.52` | `-2.85%` |
| saturated samples | 0 | 0 | no change |

At `tau ~= 15.93`, the alternate controller's exact frozen sampled poles were

```text
{0.995181633088209, 0.6712951945946475}.
```

The alternate remains locally stable and nonalternating at the selected update rate. It is modestly smoother and uses marginally less pressure, but it tracks both state components less accurately. The original index-433 eigenfunction therefore remains the better single-eigenfunction control coordinate for this target under the present cost.

This result also shows that matching an eigenvalue's autonomous frequency to target frequencies is not sufficient. A useful KLQR coordinate must encode target-versus-state error with suitable geometry and control gradient across the entire controlled trajectory. The index-433 polynomial evidently provides that geometry better despite its less obvious frequency match.

Alternate Viewer artifact:

```text
viewer_outputs/rpe-alt-krom-index431-klqr-multisine-controlled-dtauc0p001.html
```

### 17.3 Two-eigenfunction KROM experiment

The two single-function candidates were then retained together from the same original SINDy KROM, without refitting. The reduced object contains full-spectrum indices `431` and `433` and their conjugate duals, indices `432` and `434`. Thus it has two nonredundant complex eigenfunctions and four real split coordinates.

| Full-spectrum index | Continuous eigenvalue per `tau` |
|---:|---:|
| `431` | `-0.21123030699891182 + 4.822747859492991 i` |
| `433` | `-0.13029888671192555 + 5.598923003706494 i` |

The extracted KROM was immutably reserved as:

```text
archives/derived/koopman-krom-extract/
rpe-31-free-natural-period-dtau0p001-tau25-c000-t030__20260821T203315.886026Z__koopman-krom-extract__20260826T181924.655348Z__be404732.npz
SHA-256: 29a2378d6716b948b00306b4a6994a20596896b51920a6498c7fbb4897320f92
reservation label: rpe-two-eigenfunction-indices431-433
```

Everything except the KROM state dimension was held equal to the selected index-433 experiment: identical nonlinear plant and physical parameters, multisine target, equilibrium initial state, `R=1e-16`, `delta_tau=delta_tau_c=.001`, RK4, `+/-200 kPa` pressure limit, Python backend, and 25-`tau` duration. The state cost was expanded naturally from `I2` to `Q=I4`.

Controlled archive:

```text
archives/solutions/rayleigh-plesset/
rpe-equilibrium-two-eigenfunction-krom-431-433-klqr-multisine-tau25__20260826T182011.891642Z.npz
SHA-256: 1330ee074e51d73d7496c85f67e06d00bf689d8ae79484225aef361f2aec86a1
solve time: 110.777 s
```

Direct comparison with the selected single index-433 case:

| Metric | Single index 433 | Two eigenfunctions 431+433 | Two-mode change |
|---|---:|---:|---:|
| radius RMSE | `0.00635676` | `0.00513297` | `-19.25%` |
| velocity RMSE | `0.00733403` | `0.00550956` | `-24.88%` |
| max radius error | `0.01140902` | `0.00915235` | `-19.78%` |
| max velocity error | `0.02268025` | `0.01692662` | `-25.37%` |
| pressure RMS (Pa) | `23762.90` | `24062.39` | `+1.26%` |
| peak pressure (Pa) | `68792.98` | `68282.59` | `-0.74%` |
| maximum command jump (Pa) | `2919.34` | `221.06` | `-92.43%` |
| max jump, `14<=tau<=17` (Pa) | `2544.84` | `216.66` | `-91.49%` |
| pressure TV, `14<=tau<=17` (Pa) | `264171.59` | `253285.84` | `-4.12%` |
| saturated samples | 0 | 0 | no change |

At `tau ~= 15.93`, the frozen continuous closed-loop poles per `tau` were

```text
-474.984929
-5.1874243
-0.427473 +/- 5.195218 i
```

and exact-ZOH sampled poles for `delta_tau_c=.001` were

```text
0.5237163908
0.9948262850
0.9995591297 +/- 0.0051929735 i
```

All sampled poles are inside the unit circle. In the dimensionless/model-scaled coordinates at this state, the controllability matrix had singular values

```text
7.070686755e-4, 1.320591182e-4, 3.506664854e-6, 6.935121844e-7
```

and numerical condition number `1019.55`, with full numerical rank four. This is substantially worse conditioned than either two-real-coordinate single-eigenfunction realization. A controllability matrix formed directly from physical-time powers appears rank deficient and has condition number `3.78e19`; that diagnostic is unit/scaling sensitive and should not be interpreted as contradicting the full-rank dimensionless result.

The two-function controller is the best tested candidate so far: it materially improves both tracking coordinates while almost eliminating the large adjacent command jumps, at the cost of a `1.26%` increase in pressure RMS and a more poorly conditioned control realization. It still recorded 25,296 command changes and 148 same-time changes, so the floating controller-clock defect remains. The result should be retained as a research candidate and rerun after the integer-tick clock correction before it is considered a production default.

Viewer artifact:

```text
viewer_outputs/rpe-two-eigenfunction-krom-431-433-klqr-controlled.html
```

### 17.4 Off-equilibrium, large-amplitude target and reference-feedforward diagnosis

To distinguish KROM inadequacy from controller-architecture error, the two-eigenfunction KROM was tested on a harder target that remains inside the radius range of its source trajectory:

```text
r_ref(tau) = 1.5 + 0.2 sin(4 tau - pi/2) + 0.3 sin(3 tau - pi/2)
u_ref(tau) = d r_ref / d tau
```

The phases make the target start at the plant equilibrium, `(r,u)=(1,0)`, while its long-time center is `r=1.5`. Over `0<=tau<=25`, the target spans `1.0000<=r<=1.95066` and `-1.65768<=u<=1.65768`; the original small target spans only `0.86351<=r<=1.13649`. The KROM source trajectory spans `0.23287<=r<=2.3334`, so this experiment tests a different operating region without requiring radius extrapolation beyond the training trajectory.

Three 25-`tau` runs used the same nonlinear plant, equilibrium initial condition, two-function KROM, RK4 step `delta_tau=.001`, controller interval `delta_tau_c=.001`, and `+/-200 kPa` pressure limit:

1. the existing feedback-only law with `Q=I4`, `R=1e-16`;
2. the same law plus exact nonlinear reference feedforward, with `R=1e-16`;
3. feedforward plus one decade of feedback regularization, `R=1e-15`.

The production controller currently implements only

```text
p_c = -K(x) [psi(x) - psi(x_ref)]
```

At exact reference matching this necessarily returns zero control. A moving RPE reference generally requires nonzero pressure, so exact matching is not an invariant solution of the implemented closed loop. For this experiment the required reference pressure was computed by nonlinear inversion of the dimensionless RPE acceleration equation:

```text
p_ff = P0 r_ref [d u_ref/dtau - f_2(r_ref,u_ref,p_c=0)]
p_c  = p_ff - K(x) [psi(x) - psi(x_ref)]
```

This feedforward pressure has mean `81.407 kPa`, RMS `90.578 kPa`, and peak magnitude `140.912 kPa`; it never reaches the pressure limit. For comparison, the earlier small target itself requires a nonzero feedforward with RMS `25.471 kPa` and peak magnitude `71.741 kPa`. Thus its imperfect matching was structurally expected even before accounting for KROM approximation.

Measured results:

| Metric | Feedback only, `R=1e-16` | Feedforward, `R=1e-16` | Feedforward, `R=1e-15` |
|---|---:|---:|---:|
| radius RMSE | `0.01096827` | `0.000473329` | `0.00000179399` |
| velocity RMSE | `0.01239878` | `0.001426797` | `0.00000950702` |
| radius bias, `5<=tau<=25` | `-0.00971124` | `-0.000387657` | `8.12e-8` |
| max radius error | `0.01985195` | `0.00121740` | `0.00000678329` |
| max velocity error | `0.03144065` | `0.00679165` | `0.0000436584` |
| pressure RMS (Pa) | `114538.55` | `114837.86` | `90576.83` |
| peak pressure (Pa) | `200000` | `200000` | `140913.44` |
| saturated samples | `3419` | `3365` | `0` |
| maximum adjacent jump (Pa) | `400000` | `400000` | `467.94` |
| pressure total variation (Pa) | `1.7777e9` | `1.7422e9` | `3.9381e6` |
| pressure TV, `14<=tau<=17` (Pa) | `8.7834e6` | `8.6962e6` | `188797.25` |
| solve time (s) | `130.94` | `136.70` | `367.13` |

Feedforward at unchanged `R` reduced radius RMSE by `95.7%` and velocity RMSE by `88.5%`, proving that missing reference motion—not the selected eigenfunctions—was the primary source of tracking offset. It did not cure chatter because the feedback gain remained excessive. Raising `R` to `1e-15` reduced the gain norm from approximately `1.66e8` to `5.25e7`; with feedforward carrying the nominal motion, this regularization improved rather than degraded tracking. The final result is nearly exact, nonsaturating, and smooth.

Frozen sampled poles remain inside the unit circle in all three cases, so the harder-target chatter is not the earlier negative sampled-pole mechanism. Along the final trajectory, sampled spectral radius stays near `0.999573` and the fastest real sampled pole lies between `0.495` and `0.918`. The remaining risk is numerical: the four-real-coordinate controllability condition number remains about `1001` to `1036`, and the floating clock still records 25,296 command changes with 148 same-time changes.

Immutable solution archives:

```text
archives/solutions/rayleigh-plesset/
rpe-equilibrium-two-krom-431-433-klqr-offset1p5-large-tau25__20260826T182951.921418Z.npz
SHA-256: 5f4462896812b9e89279c0f0d1e3027cf45e9b4b50a08df1096dbce06f1d7d15

rpe-equilibrium-two-krom-431-433-klqr-offset1p5-large-feedforward-tau25__20260826T183619.663043Z.npz
SHA-256: 5dd7d13164ae58b41af8565c4603c98ec23c5d5a755ab31113a5900cf76891de

rpe-equilibrium-two-krom-431-433-klqr-offset1p5-large-feedforward-R1e-15-tau25__20260826T184121.490182Z.npz
SHA-256: d7d9249661ed69b3461c9882298b62dd1d4ed561ac50702e47b728344ae5a889
```

Viewer artifacts:

```text
viewer_outputs/rpe-two-krom-offset1p5-large-feedback-only.html
viewer_outputs/rpe-two-krom-offset1p5-large-feedforward.html
viewer_outputs/rpe-two-krom-offset1p5-large-feedforward-R1e-15.html
```

The feedforward controller used here was an experiment-local runtime wrapper; production controller code was not changed. The recommended implementation is a general reference-tracking contract that accepts reference state and reference derivative or a system-provided inverse-dynamics callback, combines feedforward and feedback before applying limits, and archives both components separately. It should be paired with exact discrete-time LQR/SDRE at the declared controller interval, an integer-tick sample clock, coordinate/output weighting instead of unscaled `Q=I`, and command-rate or increment penalties. KROM expansion or retraining should be considered only after those defects are corrected, or when a requested target leaves the training domain.

## 18. Recommended follow-up experiment design

Preserve the original archive and run a labeled matrix of independent cases:

1. **Isolate sampling:** keep `Q=I`, `R=1e-16`; compare controller intervals `0.005`, `0.002`, and `0.001 tau`.
2. **Isolate cost ratio:** keep `delta_tau_c=0.005`; compare `R=1e-16`, `1e-15`, and `1e-14`.
3. **Completed conservative candidate:** retain its archive as evidence that `R=1e-15` sacrifices too much tracking accuracy.
4. **Discrete design:** replace continuous CARE plus hold with an exact-ZOH discrete-time Riccati design at the declared controller interval.
5. **Normalized design:** define pressure internally as `p_c*`, use declared dimensionless state/lift coordinates, and record dimensionless `Q,R` without mixed-unit ambiguity.
6. **State-space objective:** introduce an output/state penalty if radius and radial velocity need separate semantic weights. `Q=I` in the current minimal KROM cannot distinguish them directly.
7. **Robustness:** repeat across initial target phases, target amplitudes, pressure bounds, and nearby physical parameters.
8. **Model adequacy:** compare this single-eigenfunction KROM with larger extracted KROMs. The present KROM has no state reconstruction and may omit controlled off-manifold dynamics.
9. **Clock correction:** replace the floating `floor` clock with a monotone integer-tick clock, then separately reproduce both the baseline and the selected `R=1e-16`, `delta_tau_c=.001` case before combining the clock fix with any other tuning change.

For every run, archive and compare:

- radius and velocity RMSE and maximum error;
- pressure mean, RMS, peak, and total variation;
- pressure-update jump distribution;
- saturation count and duration;
- frozen continuous and exact sampled poles along the trajectory;
- `||K||`, `||B_psi||`, and controllability conditioning;
- solve time and number of Riccati solves.

## 19. Faithful reproduction recipe

The following expresses the original API-level construction. Running it creates a new immutable solution; do not reuse the original label if archival naming policies change.

```python
from pathlib import Path
import numpy as np

from control import KoopmanLQRController, rpe_control_target
from rom import load_krom_archive
from solver import SolveRequest, SolverOptions, solve

root = Path.cwd()
krom_path = root / (
    "archives/derived/koopman-krom-extract/"
    "rpe-31-free-natural-period-dtau0p001-tau25-c000-t030__"
    "20260821T203315.886026Z__koopman-krom-extract__"
    "20260826T164339.807992Z__ebb5aedc.npz"
)

krom = load_krom_archive(krom_path)
system = krom.system
T0 = system.dimensionless_parameters.time_scale_s

controller = KoopmanLQRController(
    krom,
    Q=np.eye(2),
    R=np.array([[1e-16]]),
    target=rpe_control_target("rpe_multisine_equilibrium_small"),
    target_space="state",
    control_limit=np.array([200000.0]),
)

request = SolveRequest(
    system=system,
    kind="ode",
    time_span=(0.0, 25.0*T0),
    initial_state=np.array([system.equilibrium_radius_m, 0.0]),
    controller=controller,
    options=SolverOptions(
        method="RK4",
        rtol=1e-8,
        atol=1e-10,
        step_size=0.001*T0,
        control_sample_time=0.005*T0,
        use_dimensionless_form=None,
        execution_backend="python",
    ),
    label="rpe-equilibrium-minimal-krom-klqr-multisine-tau25-reproduction",
)

result = solve(request)
assert result.status == "complete"
print(result.archive_path)
```

Before interpreting differences from the original run:

1. verify the minimal-KROM SHA-256;
2. verify the relevant source hashes or use a preserved checkout containing those exact files;
3. verify Python, NumPy, and SciPy versions;
4. compare the new archive metadata field by field;
5. compare numerical metrics before comparing plots.

## 20. Interpretation boundaries

- The sampled-pole analysis freezes the state-dependent `B_psi` and gain at each point. It diagnoses the local mechanism but is not a theorem about the complete nonlinear, time-varying closed loop.
- The single complex eigenfunction provides two real lifted coordinates but no direct state reconstruction. Good tracking in this run does not establish global KROM validity.
- Continuous CARE stability does not guarantee stability or monotone decay after sample-and-hold implementation.
- The preserved clock is not perfectly periodic: it contains 29 duplicate same-time command changes and recurring half-step phase offsets. Reproducing the scientific configuration with a corrected clock will not be bitwise identical to the baseline.
- Identity `Q` is basis-dependent. Rescaling or rotating an eigenfunction changes the meaning of `Q=I` unless the cost is transformed consistently.
- The very small dimensional `R` and the ill-conditioned controllability matrix make numerical scaling part of the control-design problem, not merely an implementation detail.
- The original archive is preserved evidence. Follow-up tuning must produce separately named archives and reports.

## 21. Case-study lesson

This case is valuable because every individual component appears reasonable in isolation:

- the continuous KROM generator is stable;
- the frozen continuous-time LQR poles are stable;
- the pressure command stays far below saturation;
- state tracking errors are small;
- the numerical integrator uses a finer step than the controller interval.

Nevertheless, the implemented sampled controller contains a strongly alternating mode. The lesson is that KLQR/SDRE validation must include the implemented controller clock and exact held-input dynamics. Continuous Riccati poles alone are insufficient.

## 22. MATLAB `soft` parity audit and revised diagnosis

This section supersedes the cross-suite interpretation in Section 17.4. The
feedforward experiments and metrics in that section remain valid observations
about the Python controller, but the MATLAB source does **not** contain an
equivalent automatically applied reference-pressure feedforward term.

### 22.1 What MATLAB `soft` actually means

The MATLAB compute path is direct and state dependent:

1. `compute_node.compute("soft")` calls `ode15s` with `advance` as its ODE
   right-hand side (`compute_node.m`, lines 100--101).
2. Every invocation of `advance` calls `ctrl.u(x,t)` before evaluating the
   plant equation (`compute_node.m`, lines 236--244).
3. KLQR constructs the gain with
   `lqr(pl.dyn.A,obj.B(x,t,pl,o),Q,R)` (`CTRL.m`, line 470).
4. For a KROM plant, `obj.B(x,t,pl,o)` evaluates the eigenfunction gradient,
   the observable Jacobian, and the current physical input map
   (`plant.m`, line 79).
5. The tracking command is only `u=-K*K_err` (`CTRL.m`, lines 137 and 146).

Consequently, MATLAB recomputes the current physical input matrix, lifted
`B_psi`, continuous-time CARE/LQR gain, and control command at every internal
`ode15s` right-hand-side evaluation. This includes trial states and evaluations
that may later be rejected by the adaptive solver. There is no independent
controller clock and no zero-order command hold in the `soft` path. The saved
time spacing is output sampling, not the feedback-update interval.

MATLAB's `hard` path is different: its custom RK4 implementation computes one
command before the four RK stages and holds that command through the stages.
It is therefore closer to a fixed-step sampled controller than `soft` is.

### 22.2 Difference from the current Python sampled implementation

The ordinary Python solve recomputes `B_psi` and the CARE solution when the
controller clock requests a new command. Between those ticks, it holds the
scalar pressure command. The nonlinear RPE input map is still reevaluated at
each Runge--Kutta stage while applying that held command, but the lifted input
map and gain are not. The distinction is:

```text
MATLAB soft:    RHS state -> B(x) -> B_psi(x) -> CARE -> u -> RHS
Python sampled: controller tick -> B_psi(x) -> CARE -> held u
                every RHS stage -> B(x_stage) * held u
```

This is large enough to change behavior for an SDRE/KLQR controller. Earlier
experiments already showed sensitivity to a nominal `0.001 tau` versus
`0.005 tau` controller period. MATLAB `soft` goes further: its effective update
times are the adaptive solver's internal evaluations, not a regular sampled
period. The Python floating-point clock defect, which can produce duplicate
same-time controller updates, is another independent mismatch and must not be
mistaken for continuous state dependence.

### 22.3 No hidden pressure feedforward

The source audit found no automatic `p_ff`, inverse-dynamics pressure, target
derivative, or integral tracking term in the MATLAB KLQR command. The quantity
named `FF` in the Bubble module is a far-field acoustic-echo observable, not a
controller feedforward contribution. Therefore:

- Python feedforward substantially improved the tested Python trajectory;
- it is a useful two-degree-of-freedom control design candidate;
- it does not explain the reported MATLAB tracking result; and
- it must remain disabled in a strict parity experiment.

The prior claim that missing feedforward was the primary cross-suite cause was
too strong. The more plausible parity differences are feedback evaluation
policy, KROM selection/size and normalization, target lifting, tuning/scaling,
and the plant input sign convention.

### 22.4 Plant input sign convention

MATLAB Bubble defines its dimensionless pressure input with a negative sign,
`B(r)=-[0;Eu/r]` (`Bubble.m`, line 167): positive external acoustic pressure
compresses the bubble. The current Python RPE input map uses the opposite sign.
This must be reconciled before comparing archived pressure histories.

If both the controller model and plant input matrix are sign-flipped
consistently, ideal pure-LQR closed-loop state dynamics are invariant: the gain
and command reverse sign together, so `B*u` is unchanged. Thus the mismatch is
important for physical command semantics, feedforward construction, limits,
and archive interpretation, but by itself it does not explain the tracking
gap of a consistently modeled feedback-only controller.

### 22.5 KROM mismatch is likely material

The active MATLAB examples do not establish that the controller used the same
two eigenfunctions selected for the Python experiment. The MATLAB workflow
creates a full H-DMD model with a `10000 x 2000` delay matrix, optionally fits
SINDy eigenfunctions, and then applies observable-based extraction with
thresholds and caps (for example, `threshold=[0.99,0.99,0.99]` and
`caps=[30,3,200]`). The control scripts refer to workspace variable `kROM{1}`;
several two-mode extraction statements are commented out. A claim of MATLAB
parity therefore requires preserving the exact MATLAB KROM membership,
normalization, observable map, spectrum convention, and tuning—not merely a
two-function Python KROM with nearby eigenvalues.

This also explains an otherwise striking result from the feedback-only trial
against the MATLAB scanner's fixed target `[2,0]`: the selected Python
two-function KROM generated an almost constant small command and left the
bubble near equilibrium. That lifted error/control direction is poorly suited
to the requested state displacement; it is not evidence that continuous SDRE
is intrinsically ineffective.

### 22.6 Feedback-only parity experiments attempted

All trials below started the controlled RPE at equilibrium, used no pressure
feedforward, and did not mutate MATLAB or production Python source.

- The active hard MATLAB waveform (offset `4`, amplitudes `[1.5,1]`, angular
  frequencies `[4,3]`) was attempted with continuous per-RHS or per-RK-stage
  gain recomputation. Aggressive MATLAB-like and relaxed `R` values produced
  enormous initial commands or drove the model outside positive-radius/
  finite-state validity before completing the requested horizon. Adaptive
  stiff trials were also prohibitively slow because each internal evaluation
  solved a CARE problem.
- The commented MATLAB single-sine case (offset `2`, amplitude `1.5`, angular
  frequency `4.59333`) likewise left the valid domain under the two-function
  KROM.
- The explicit scanner fixed target `[2,0]`, with continuous four-stage RK4
  recomputation, `Q=I`, `R=1e-8`, `delta_tau=0.005`, and horizon `10 tau`, did
  complete. It made 8,000 controller evaluations in 175.24 seconds but failed
  to move the radius meaningfully: radius RMSE was `1.000039`, final radius was
  `0.9999817`, pressure RMS was `16.973 Pa`, and peak pressure was `16.988 Pa`.

Completed immutable archive:

```text
archives/solutions/rayleigh-plesset/
rpe-continuous-stage-fixed-r2-feedback-only-Q1-R1e-8-tau10__20260826T194822.218528Z.npz
SHA-256: 1bb49733d2679b3c790015468cb01b95714c1c598420606a5a8982fe202533b8
```

The failed waveform experiments created no accepted solution archive.

### 22.7 Feasibility and performance of a closer formulation

A closer MATLAB-soft solve is feasible, but a literal Python implementation is
too expensive for routine work. The main cost is not evaluating `B(x)` or
eigenfunction gradients; it is repeatedly solving the continuous CARE during
an adaptive stiff integration. The unsuccessful trials took minutes while
advancing little physical time, and one fixed-target run required about 175
seconds for only `10 tau`.

The parity implementation should therefore be a named, testable policy rather
than an implicit change to the default solver:

1. Add `control_evaluation="continuous_rhs"` alongside the ordinary sampled
   policy. It must be stateless at trial evaluations and must not write every
   rejected solver call into the scientific control history.
2. Record accepted/output-time control history separately from an optional
   diagnostic evaluation tape.
3. First reproduce MATLAB with the exact KROM, target, sign, dimensionless
   `Q/R`, tolerances, and initial state.
4. Then accelerate CARE evaluation with warm starts, state-grid interpolation,
   gain scheduling, or a structure-aware low-dimensional solver. Cache only
   when a declared state tolerance makes the approximation explicit.
5. Compare the accelerated policy against literal continuous-RHS recomputation
   on short horizons before using it for scientific runs.

### 22.8 Startup pressure spikes: accepted for parity, deferred for design

Large initial pressure spikes are expected when the bubble begins at
equilibrium while the target begins far away or with a large nonzero velocity,
especially with an extremely small control penalty. They were also observed
and disliked in the MATLAB workflow. They are permitted temporarily so strict
parity work does not silently change the controller problem.

Once parity is established, spike suppression must be studied independently.
Candidate mechanisms are a smooth reference/engagement ramp, bumpless transfer,
command magnitude and slew-rate constraints, control-increment penalties,
anti-windup if integral action is added, and finite-horizon preview. Each must
be evaluated for its effect on post-transient tracking; clipping the first
sample without analysis is not an acceptable fix. Preserve both the parity
baseline and every spike-mitigated run as separate archives.

## 23. Revised immediate recommendation

Do not tune `Q/R` further against the current two-function controller and infer
MATLAB parity from the result. First reconstruct the exact MATLAB-controlled
KROM and implement a diagnostic continuous-RHS evaluation policy with correct
pressure semantics. Then compare, in order:

1. exact MATLAB-like continuous feedback with no feedforward;
2. regular sampled feedback using the same KROM and cost;
3. accelerated continuous/gain-scheduled feedback;
4. optional feedforward as a deliberate Python extension; and
5. startup-spike mitigation after the parity baseline is stable.

This ordering separates implementation parity, model adequacy, numerical
acceleration, and control-design improvement instead of allowing one to mask
another.

## 24. Implemented solver policy boundary

The solver now exposes `SolverOptions.control_evaluation` with three declared
values:

- `"sampled"` requires `control_sample_time` and uses an integer-tick
  zero-order hold;
- `"continuous_rhs"` forbids a sample time and recomputes an RHS-safe
  controller at every solver evaluation; and
- `"automatic"` preserves compatibility by selecting between those policies
  from the presence of `control_sample_time`.

Continuous-RHS mode rejects stateful controllers and estimators. It archives
the complete internal controller-evaluation tape separately from the control
series sampled on accepted solution states. Sampled mode now records only real
controller updates rather than recording every reuse of a held value. The
integer-tick clock also removes the known `floor(t/step)` duplicate-update
mechanism. This is an implementation capability, not yet evidence of MATLAB
parity; the exact KROM, input convention, target, and tuning still have to be
matched in the bounded diagnostic run.

## 25. Bounded continuous-versus-sampled diagnostic

The reproducible command module `scripts.rpe_klqr_parity_diagnostic` was run
for only `0.05 tau`, with `delta_tau=0.001`, the current two-function KROM,
the offset-4 mixed target, `Q=I4`, `R=1e-8`, equilibrium initial state, no
feedforward, and no control limit. This short horizon intentionally retains the
startup spike while avoiding an invalid long extrapolation.

| quantity | continuous RHS | sampled, `delta_tau_c=0.001` |
|---|---:|---:|
| controller evaluations | 200 | 51 |
| initial pressure | `-66,157.604 Pa` | `-66,157.604 Pa` |
| peak absolute pressure | `187,761.340 Pa` | `186,915.280 Pa` |
| final radius | `0.988931003` | `0.989029297` |
| final radial velocity | `-0.534280815` | `-0.529118568` |

The maximum difference between solution-aligned pressure histories was
`846.061 Pa`. The continuous formulation compressed the bubble slightly more
over this initial interval. Both policies initially move away from the target,
which starts at radius `4` with positive radial velocity. Therefore continuous
state dependence is real and measurable, but it does not repair the selected
two-function KROM's unfavorable initial lifted control direction. This further
supports reconstructing the actual larger/observable-selected KROM before
interpreting `Q/R` tuning or long-horizon parity.

Immutable archives:

```text
rpe-continuous-rhs-feedback-only-rpe_multisine_offset_4_mixed-Q1-R1e-08-tau0-05__20260826T200350.791582Z.npz
SHA-256: 3a8ccca5c9e2fc7d69203cbc7f0774ba20654a41cb5be8601714335dba8ee9b4

rpe-sampled-feedback-only-rpe_multisine_offset_4_mixed-Q1-R1e-08-tau0-05__20260826T200416.124984Z.npz
SHA-256: 3ce955f9400d275bb37829f50e39462cf59449c14fbc22bb51efc7a283bb5d2e
```

## 26. Versioned correction of RPE pressure semantics

The Python RPE input has now been corrected and versioned across dimensional,
dimensionless, analytic-Jacobian, KROM, and native-adapter paths:

```text
external-acoustic-pressure/v2:       positive pressure compresses the bubble
legacy-positive-acceleration/v1:     positive input expands the bubble
```

Newly constructed systems use the physical `v2` convention. An archive that
lacks the field is reconstructed as legacy `v1`, preserving all earlier
solution and KROM behavior. Existing archives were not rewritten. An analytic
KROM may be rebound to a control-equivalent revised system only through the
explicit `with_analytic_control_system(..., reason=...)` operation. That
operation copies the KROM, preserves its autonomous spectral source, rejects
fitted linear/bilinear control matrices, validates system type/state
coordinates, and records both system snapshots and the reason.

The bounded diagnostic was repeated with that explicit physical rebinding. As
predicted for a consistent pure LQR sign transformation, the state trajectories
and pressure magnitudes were unchanged to reported precision, while every
pressure command reversed sign:

| quantity | physical continuous RHS | physical sampled |
|---|---:|---:|
| initial external pressure | `+66,157.604 Pa` | `+66,157.604 Pa` |
| peak absolute pressure | `187,761.340 Pa` | `186,915.280 Pa` |
| final radius | `0.988931003` | `0.989029297` |
| final radial velocity | `-0.534280815` | `-0.529118568` |

This resolves the physical interpretation bug but does not resolve target
tracking: the selected two-function KROM still commands compression at startup.
The result proves that a consistently applied input-sign correction cannot by
itself change pure-feedback closed-loop state dynamics. The remaining leading
hypothesis is KROM membership/observable selection, followed by scaling and
cost definition—not hidden feedforward or pressure-sign inconsistency.

Corrected-convention immutable archives:

```text
rpe-continuous-rhs-feedback-only-rpe_multisine_offset_4_mixed-Q1-R1e-08-tau0-05__20260826T201215.259694Z.npz
SHA-256: e346d2e88428142469a9e6a023f29f218f698f7cbc56849e1eb9c584103c32f5

rpe-sampled-feedback-only-rpe_multisine_offset_4_mixed-Q1-R1e-08-tau0-05__20260826T201218.765813Z.npz
SHA-256: adafa418403eedce7dd0b3fdeab73b4e55196cc0f2100e288c0a5e8766ff07bc
```

## 27. Observable-mode selection and validity-domain diagnosis

A grouped state-observability selector was applied to all 1,001 nonredundant
modes of the 2,000-entry SINDy spectrum. Because every polynomial
eigenfunction is a linear combination of only 15 dictionary features, the
selector operates on feature Gram matrices rather than loading the archived
`25001 x 2000` complex eigenfunction array. It completed in approximately
three seconds and selected five complex modes before meaningful independent
improvement was exhausted.

Selected full-spectrum indices and continuous eigenvalues:

| order | full index | continuous eigenvalue |
|---:|---:|---:|
| 1 | 433 | `-0.130299 + 5.598923 i` |
| 2 | 266 | `-0.438168 + 764.490786 i` |
| 3 | 767 | `-7.678146 + 953.067918 i` |
| 4 | 1983 | `-6.899675 + 3141.592654 i` |
| 5 | 272 | `-0.709514 + 809.443192 i` |

Mode 433—the previously preferred single eigenfunction—was selected first and
had the greatest equilibrium control authority. This independently validates
that earlier minimal-mode choice. However, the final linear state-observation
scores were only:

```text
radius R^2:          0.3454437181
radial-velocity R^2: 0.8154874705
stop reason:         minimum relative improvement
```

Thus the present SINDy eigenfunction family does not provide a high-fidelity
linear state-observation map, even when more modes are retained. It cannot
reproduce the MATLAB observable-extraction result merely by changing the
number of retained modes.

The validity-domain audit found a more fundamental mismatch:

```text
training radius:       0.23287463 <= r <= 2.3334
training velocity:   -26.90988407 <= u <= 26.91080092
hard target at tau=0:  r=4, u=9
radius excursion:      79.3421% of the complete training-radius span above max
```

The hard offset-4 target is already outside the polynomial training domain at
startup and moves farther outward. All such runs are extrapolation studies,
not KROM validation. KROM archives now derive and expose their training domain;
KLQR/KMPC default to an error outside it. The parity diagnostic used the
explicit `target_domain_policy="ignore"` override so the failure mechanism
could still be reproduced and labeled.

Over the same physical-sign, continuous-RHS, feedback-only `0.05 tau` test:

| quantity | prior two-function KROM | observable-selected five-mode KROM |
|---|---:|---:|
| initial pressure | `66,157.604 Pa` | `15,048.940 Pa` |
| peak absolute pressure | `187,761.340 Pa` | `37,900.144 Pa` |
| final radius | `0.988931003` | `0.997603715` |
| final velocity | `-0.534280815` | `-0.112373602` |

The five-mode KROM reduces the erroneous compression substantially but does not
reverse it. Adding observable/control-authority modes therefore improves the
startup transient without solving tracking. The next scientifically valid
step is to train a KROM whose free-response domain covers the intended target
radius range and whose extraction includes a validated state-observation map.
`Q/R` tuning against the present offset-4 extrapolation is not meaningful.

Immutable selected KROM:

```text
rpe-31-free-natural-period-dtau0p001-tau25-c000-t030__20260821T203315.886026Z__koopman-krom-extract__20260826T202020.762433Z__bdca807e.npz
SHA-256: e32d873bbc3a3c38bdd6f2c89df877938cee424dea57beadb57c685ed33a45d4
```

Explicit-extrapolation diagnostic archives:

```text
two-function: rpe-continuous-rhs-feedback-only-rpe_multisine_offset_4_mixed-Q1-R1e-08-tau0-05__20260826T202241.495608Z.npz
SHA-256: fc8e299c6871408908b167af3560de17c21b0fc69974aced76dc6987ca9f1c7d

five-mode: rpe-continuous-rhs-feedback-only-rpe_multisine_offset_4_mixed-Q1-R1e-08-tau0-05__20260826T202244.811659Z.npz
SHA-256: 3d432285c53291dd35f261756ab3c40830ae58c0d537d22bd3d6923791124026
```

## 28. Startup-spike parity policy

Large initial commands are retained while continuous state-dependent KLQR
parity is being established. The diagnostic baseline uses
`control_limit=None`, with no clipping, command ramp, slew limiter, or target
blending. This intentionally preserves the physically undesirable startup
spike so controller formulation, KROM membership, and state-dependent gain
updates can be compared without a hidden transient treatment.

Spike suppression is a deferred and separately versioned control-design study.
Candidate mechanisms include bumpless transfer, reference governors, smooth
target engagement, magnitude and slew constraints, control-increment cost, and
preview control. A mitigation must reduce peak pressure and pressure slew
without degrading post-engagement tracking error, phase accuracy, or stability
margins. Its results must be archived separately from the parity baseline.

## 29. State-observation validation and next training domain

The five-mode observable-selected KROM was fitted with an explicit affine map
from its ten real/imaginary lifted coordinates to `(r,u)`. A deterministic 20%
interleaved holdout was excluded from fitting. The map failed promotion:

| state | fit R2 | held-out R2 | held-out normalized RMSE |
|---|---:|---:|---:|
| radius | `0.345551` | `0.344244` | `0.169835` |
| velocity | `0.815555` | `0.814929` | `0.027583` |

Its normalized condition number was `1.3244e16`, above the `1e12` gate. No
observation-map KROM archive was published. This replaces the earlier informal
state-observability score with an executable, held-out capability contract.

The canonical RPE target envelopes require dimensionless radius through
`6.833333` and velocity through `+/-19.373155`. The existing KROM reaches only
radius `2.3334`, though its velocity coverage is already wider. The dry-run
coverage planner therefore proposes a radius-only expansion with 10% outward
margin to `7.493379`, preserves the existing lower/velocity bounds, and emits
72 new candidate initial conditions after removing points already covered.
Nothing was scheduled. The plan is gated for screening because extreme free
responses may approach the RPE `r=0` singularity.

## 30. Retained offset-1.5 clean-pathway follow-up

The later clean-isolation rerun and adjacent-eigenfunction comparison are
documented separately in
`RPE_OFFSET1P5_MINIMAL_KROM_FAILURE_CASE_STUDY_2026-08-27.md`. That report is
the authoritative record for the unchanged-target, continuous-RHS,
`Q=5000 I2`, `R=1e-4`, uncapped feedback-only failure. It identifies the
hidden initial acceleration demand, severe target-path eigenfunction-map
folding, physical-to-lifted error compression, and incorrect control-pressure
phase while ruling out gross training-box extrapolation, loss of `B_psi`,
lifted uncontrollability, CARE failure, and insufficient command RMS.
