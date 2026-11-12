# Colbrook/coauthor native adaptation provenance

All 14 publicly accessible reference repositories in the
[pinned manifest](../../tools/experiments/miscellaneous_computes/references/colbrook/repositories.json)
were downloaded on 2026-11-11 to ignored `build/reference_sources/colbrook`.
External datasets and submodules were not downloaded. No upstream installers
were executed. Selected mpEDMD and Rigged DMD routines were executed only on the
synthetic validation fixture.

The native implementation uses DynSysLab's compiled numerical backend and owns
its data, serialization and orchestration contracts. The following mathematical
routines were consulted/adapted; preserve these upstream license notices in
source and binary distributions:

- ResDMD residual and minimum-residual formulations:
  [BSD-2-Clause notice](colbrook/MColbrook--Residual-Dynamic-Mode-Decomposition--0f7f1bf67832-LICENSE.txt).
  Native files: `KoopmanSpectralFit.cpp`, `KoopmanResiduals.cpp`, `KoopmanStochastic.cpp`.
- mpEDMD Procrustes (`mpEDMDqr.m`):
  [BSD-2-Clause notice](colbrook/MColbrook--Measure-preserving-Extended-Dynamic-Mode-Decomposition--e473dee7b32d-LICENSE.txt).
  Native file: `KoopmanSpectralFit.cpp`; SVD whitening replaces the reference QR basis.
- Rigged DMD and rational smoothing (`riggedDMD.m`, `rational_kernel.m`):
  [BSD-2-Clause notice](colbrook/MColbrook--Rigged-Dynamic-Mode-Decomposition--013f5a429868-LICENSE.txt).
  Native file: `KoopmanMeasures.cpp`; direct resolvent solves avoid depending on an
  eigendecomposition for the smoothing calculation.
- MultDMD indicator objective (`construct_basis.m`):
  [MIT notice](colbrook/NBoulle--MultDMD--85277b046a9a-LICENSE.txt).
  Native file: `KoopmanSpectralFit.cpp`; deterministic farthest-point centers
  replace the example's k-means dependency.
- SpecRKHS Gram/residual construction (`generate_matrices.m`, `verified_eigenvalues.m`):
  [MIT notice](colbrook/GustavConradie1--SpecRKHS--8b888786c3e7-LICENSE.txt).
  Native file: `KoopmanRkhs.cpp`; real Gaussian kernels currently supported.

Principal angles, Fejer autocorrelation smoothing, weighted regression, TLS,
dictionary construction, archive/reference contracts and personalized PageRank
are original implementations of the published mathematical constructions.
No unlicensed TrustKoopman, adversarial-example or ENSO source was copied.
No MATLAB source is embedded in the runtime, and MATLAB/Python is not required
to execute the native methods. The MATLAB parity script is optional development
evidence. These implementation statements do not claim all theorems or complete
workflows from the cited repositories have been ported.

The [architecture record](../architecture/SPECTRAL_KROM.md) lists implemented
calculations, assumptions, tests, and remaining extensions explicitly.
