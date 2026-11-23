"""Periodic-family spectral comparison through the compiled post-processing APIs.

This layer orchestrates registered analytic sampling; it contains no Fourier,
Koopman fitting, quadrature, filtering or error algorithms. Numerical results
belong in caller-selected archives/research folders, not this package.
"""
from __future__ import annotations

from dynsyslab_client.periodic import sample, spectral_sweep, spectral_compare, quadrature_audit, fit_delay_model, audit_delay_forecast, packet_features


def sample_bank(system, initial_states, observables, *, weights, samples,
                maximum_harmonic, executable=None, timeout=120, families=None, parameters=None):
    """Sample registered periodic analytic orbits with explicitly supplied masses.

    The native sampler rejects systems/orbits lacking the requested periodic
    analytic solution. It checks their actual physical clock. This does not
    infer an invariant transverse measure or certify omitted Fourier tails.
    """
    initial_states=list(initial_states);weights=list(weights);observables=list(observables)
    if families is not None:families=list(families)
    if parameters is not None:parameters=list(parameters)
    if len(weights)!=len(initial_states):
        raise ValueError("Supply one explicit transverse weight per orbit")
    if families is not None and len(families)!=len(initial_states):
        raise ValueError("Supply one family label per orbit")
    if parameters is not None and len(parameters)!=len(initial_states):raise ValueError("Supply one parameter per orbit")
    bank=[]
    for i,(initial,weight) in enumerate(zip(initial_states,weights)):
        request=dict(system=system,initial_state=list(initial),observables=list(observables),
                     samples=samples,maximum_harmonic=maximum_harmonic)
        orbit=sample(request,executable=executable,timeout=timeout)
        bank.append(dict(channels=orbit["fourier_series"],weight=weight,
            parameter=parameters[i] if parameters is not None else None,
            family=families[i] if families is not None else "all",initial_state=list(initial),
            source=dict(request=request,equations_fingerprint=orbit.get("equations_fingerprint"),
                        period=orbit["period"],resolved_parameters=orbit.get("resolved_parameters"),
                        dimensionless_groups=orbit.get("dimensionless_groups"),
                        coefficient_source="Native FFT quadrature of registered analytic physical-time orbit")))
    quadrature_audit(dict(training=bank,audit_phase_counts=[]),executable=executable,timeout=timeout)
    return bank


__all__=["sample_bank","spectral_sweep","spectral_compare","quadrature_audit","fit_delay_model","audit_delay_forecast","packet_features"]
