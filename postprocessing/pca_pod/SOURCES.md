# PCA/POD provenance and optional external suites

This layer vendors no third-party source. The four algorithms are compact local
implementations built from SVD and principal-component-pursuit formulations.

Reviewed optional external references (all remain optional rather than merged):

1. [scikit-learn IncrementalPCA](https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.IncrementalPCA.html)
2. [pyMOR](https://pymor.org/) — POD and model-reduction framework
3. [PySPOD](https://pyspod.readthedocs.io/en/latest/) — spectral POD
4. [snapshot POD reference implementation](https://github.com/vassilikitsios/snapshot_pod_rom_py)
5. [pyLOM](https://arnaumiro.github.io/pyLowOrder/index.html) — HPC POD/DMD

For very large data, install a chosen external backend separately rather than
copying its code into dynsyslab.
