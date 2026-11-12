# Post-processing dependencies

Post-processing remains independent of `solver/` and `dynsyslib/` by default.

The sole exception is `koopman/koopman_sindy.py` in **domain mode**. That branch
performs a lazy import of `dynsyslib.interfaces.ODESystem` solely to evaluate the
provided dynsyslib system's equations of motion on user-specified domain samples.
Trajectory-based Koopman-with-SINDy analysis does not import `dynsyslib`.
