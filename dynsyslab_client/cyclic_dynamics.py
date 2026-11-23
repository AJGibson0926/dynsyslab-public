"""Thin orchestration of native local state-pair fits over a finite cycle."""
from __future__ import annotations


def fit_cyclic_state_pair_models(workflow, phase_pairs, *, sample_interval=1.0,
                                 svd_tolerance=1e-10, validation_pairs=None):
    """Fit one unconstrained linear DMD model per known cycle phase.

    Supply centered physical deviations in sample-major native data records
    (dimension, samples, x, y). The caller supplies and verifies the base cycle,
    centers, phase assignments and state units. Fits use the existing native
    ResDMD API; full-rank linear features, no constant and no normalization.
    Native coordinate rollouts extract physical matrix columns, avoiding any
    assumption about internal Koopman-basis orientation or normalization.

    The output values are phase-major, row-major matrices for cycle_cocycle.
    This is a finite-data local approximation, not a certified derivative.
    """
    phase_pairs = list(phase_pairs)
    if not phase_pairs:
        raise ValueError('At least one phase is required')
    dimension = phase_pairs[0]['dimension']
    if not isinstance(dimension, int) or not 1 <= dimension <= 16:
        raise ValueError('Physical dimension must be an integer from 1 to 16')
    if any(pair['dimension'] != dimension for pair in phase_pairs):
        raise ValueError('Every phase must have the same physical dimension')
    if validation_pairs is not None and len(validation_pairs) != len(phase_pairs):
        raise ValueError('Supply one validation record per phase')
    matrices, models, validation = [], [], []
    for phase, data in enumerate(phase_pairs):
        fit = workflow.koopman_spectral(dict(action='fit', data=data, settings=dict(
            method='resdmd', rank=dimension, svd_tolerance=svd_tolerance,
            sample_interval=sample_interval,
            dictionary=dict(kind='linear', constant=False, normalize=False))))
        if fit['rank'] != dimension:
            raise ValueError(f'Phase {phase} has insufficient independent training directions')
        model = fit['model']; columns = []
        for component in range(dimension):
            initial = [float(j == component) for j in range(dimension)]
            rollout = workflow.koopman_spectral(dict(action='rollout', model=model,
                reconstruction='coordinates', initial=initial, steps=1))
            columns.append(rollout['values'][dimension:2 * dimension])
        matrices.extend(columns[j][i] for i in range(dimension) for j in range(dimension))
        models.append(model)
        if validation_pairs is not None:
            validation.append(workflow.koopman_spectral(dict(action='residuals', model=model,
                                                          data=validation_pairs[phase])))
    return dict(dimension=dimension, period=len(phase_pairs), values=matrices,
                models=models, validation=validation,
                scope='Local linear state-pair fits conditioned on supplied finite-cycle phase; no global closure or derivative certificate.')
