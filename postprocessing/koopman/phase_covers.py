"""Finite phase-cover geometry using registered native analytic flows.

This adapter composes existing queries. Exact branch arithmetic is owned by
native mathematics; trajectory evaluation is owned by the system library.
"""
from dynsyslab_client.mathematics import query
from dynsyslab_client.periodic import sample


def sample_phase_cover(system, initial_states, *, prime, samples=1024,
                       relative_frequency_tolerance=1e-6, executable=None):
    """Sample a caller-selected prime-power orbit tower and check its clocks.

    Initial states define phase zero on each orbit. Returns uniform-time state
    samples, exact frequency-ratio labels and equation provenance. This does not
    search for orbits or prove an infinite tower exists. A failed clock check
    raises instead of silently drawing a non-equivariant covering.
    """
    states=list(initial_states)
    if not states:raise ValueError('Supply at least one orbit')
    if not 0<=relative_frequency_tolerance<=.01:raise ValueError('Invalid frequency tolerance')
    levels=query(dict(operation='phase_tower',prime=prime,digits=[0]*(len(states)-1)),executable=executable)['result']['levels']
    orbits=[];base=None
    for state,level in zip(states,levels):
        result=sample(dict(system=system,initial_state=list(state),observables=[],samples=samples),executable=executable)
        omega=result['angular_frequency']
        if base is None:base=omega
        error=abs(omega*int(level['modulus'])/base-1)
        if error>relative_frequency_tolerance:raise ValueError(f'Orbit clocks do not form the requested cover: relative error {error}')
        orbits.append(dict(level,angular_frequency=omega,period=result['period'],relative_clock_error=error,
            state_dimension=result['state_dimension'],values=result['values'],initial_state=list(state),
            equations_fingerprint=result['equations_fingerprint'],closure_error=result['closure_error']))
    return dict(schema='dynsyslab.phase-cover-geometry/v1',prime=prime,samples=samples,orbits=orbits,
        phase_convention='Uniform physical time from each supplied initial state; endpoint omitted',
        scope='Finite sampled circle coverings; not a physical-state conjugacy or an infinite-limit certificate')
