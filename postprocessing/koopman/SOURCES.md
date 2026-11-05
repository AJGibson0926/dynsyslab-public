# Koopman/DMD provenance and optional backends

This package vendors **no external implementation**. Its small, independent
implementations of DMD, Hankel DMD, and EDMD were written and tested here from
their standard linear-algebra formulations. This avoids mixing incompatible APIs
or copying source code under differing licenses.

Reviewed optional algorithm sources (well below the 20-source cap):

1. [PyDMD](https://github.com/PyDMD/PyDMD) — MIT; DMD variants including Hankel DMD.
2. [PyKoopman](https://github.com/dynamicslab/pykoopman) — Koopman/DMD framework.
3. [PySINDy](https://github.com/dynamicslab/pysindy) — sparse governing-equation discovery, complementary to Koopman analysis.
4. [datafold](https://github.com/datafold-dev/datafold) — EDMD and manifold time-series tools.
5. [pykoop EDMD documentation](https://pykoop.readthedocs.io/en/stable/pykoop.html) — EDMD pipeline reference.

External projects can be installed later as optional backends through
`adapters.py`; they are deliberately not recombined with this package's code.
