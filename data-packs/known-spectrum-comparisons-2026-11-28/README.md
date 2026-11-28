# Known Koopman spectra against mpEDMD and Rigged DMD

This pack regenerates, from nothing, the comparison of a system's exactly known continuous Koopman
spectrum with the spectrum that measure-preserving EDMD and Rigged DMD estimate from the same orbits.
It covers every catalogue system with a shell band: a phase e^{iφ} with frequency Ω(H), whose
generalized eigenfunctions δ(Ω(H(x)) − p) e^{iφ(x)} span a band of the continuous spectrum
(Mezić, *J. Nonlinear Sci.* 30, 2020, §11).

Each system runs the same job group, `jobs/10-comparison/comparison.job.json`, with its own
parameters (`jobs/00-reproduce-all.job.json`):

1. **Solve** 33 periodic orbits started at rest across the family (training) and 8 interleaved orbits
   (held out), without damping or forcing.
2. **Reference** (`band_projection`): along each training orbit, the observable's Fourier coefficients
   in the catalogue band's exact phase, by weighted Birkhoff averages (Das & Yorke, *Nonlinearity* 31,
   2018). Weighted across orbits by the canonical area dE/Ω, they give the family's spectral measure,
   smoothed by Rigged DMD's kernel (order 6, ε = 0.05), and its wave packets.
3. **Fit** (`mpedmd`, `"pool": true`): one measure-preserving EDMD fit (Colbrook, *SIAM J. Numer.
   Anal.* 61, 2023) on delays of the first state coordinate, pooled over the training orbits, with the
   same canonical-area weights and a pseudospectrum lattice.
4. **Rigged DMD** (`rigged_measure`; Colbrook, Drysdale & Horning, arXiv:2405.00782): the observable's
   density on 1024 angles and its wave packet at angle 0.5.
5. **Both on one spectrum** (`band_projection` with the krom): the reference band beside Rigged DMD's,
   with the density's relative L2 error.
6. **Held-out packets** (`packet_validation`): the fitted packets along the held-out orbits against
   the reference.
7. **Figures** (`viewer_figure`): the spectral measures (atoms on the unit circle with both
   densities), the eigenvalues over the pseudospectrum, and the packet amplitudes by orbit frequency.

## Results it reproduces

| System | Observable | Delays | Density error | Fitted / reference mass | Held-out packet error |
|---|---|---:|---:|---:|---:|
| Pendulum (librations) | angle | 32 | 9.4 % | 1.440 / 1.566 | 98 % |
| Duffing (outer orbits) | position | 32 | 25.2 % | 1.849 / 2.106 | 68 % |
| Morse | position | 32 | 24.9 % | 0.2571 / 0.2642 | 93 % |
| Pöschl–Teller | position | 32 | 14.3 % | 0.4945 / 0.4949 | 99 % |
| Cosh | position | 32 | 5.1 % | 1.689 / 1.777 | 1000 % |
| Radial twist | position | 32 | 12.6 % | 0.5931 / 0.6314 | 5.2 % |
| Kepler | q₁ | 32 | 41.8 % | 1.882 / 1.924 | 46 % |
| RPE | radius | 128 | 0.49 % | 3.710 / 3.729 | 29 % |

Density error is the relative L2 error on 1024 angles; packet error the relative L2 error of the RMS
packet amplitude at angle 0.5 over the held-out orbits. The packet centre 0.5 is the one used for the
RPE study, where the fundamental is near Ω Δt = 0.5. For most other families it falls on a weak
harmonic, where the reference packet is small and its relative error large (cosh). The RPE radius
mass, 3.710, agrees with the earlier private study's 3.7088 for the same observable. These are fits
of a fixed, modest size; the density error falls with the number of delays, as the RPE row at 128
delays shows.

## Run it

See [docs/DATA_PACKS.md](../../docs/DATA_PACKS.md). It runs in about 30 minutes on one machine and
needs a few GB of free space while it runs.
