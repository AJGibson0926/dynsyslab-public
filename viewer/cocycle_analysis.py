"""Display supplied return-map, matrix-cocycle and finite-partition results.

No integration, differentiation, spectral fitting or certification is performed
here. Native phase-space and mathematics queries own those computations.
"""
from __future__ import annotations

RED = '#cf3028'
BLUE = '#0072b2'
GOLD = '#b88400'


def multiplier_branch_figure(branch, local, *, parameter_label='Parameter'):
    """Plot preselected real multipliers; local is a separately refined branch.

    Each record provides parameter and multiplier arrays. A candidate crossing
    is displayed as a numerical bracket, never as a certified bifurcation.
    """
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4))
    for ax, record, title in zip(axes, [branch, local],
                                ['Physical periodic-orbit branch', 'Refined crossing bracket']):
        ax.plot(record['parameter'], record['multiplier'], color=BLUE, lw=2)
        ax.axhline(-1, color=RED, ls='--', lw=1.4, label='Flip candidate: multiplier = −1')
        if record.get('bracket'):
            ax.axvspan(*record['bracket'], color=GOLD, alpha=.23, label='Computed bracket')
        ax.set(xlabel=parameter_label, ylabel='Real transverse multiplier', title=title)
        ax.grid(alpha=.2)
        ax.ticklabel_format(axis='x', style='plain', useOffset=False)
        ax.tick_params(axis='x', rotation=20)
        ax.legend(fontsize=9)
    fig.suptitle('Return-map stability from the original vector field', fontsize=14)
    fig.tight_layout()
    return fig


def perturbation_figure(records, *, coordinates=(0, 1)):
    """Compare precomputed free perturbation rollouts and deviation errors.

    Records have label, time, reference_deviation and models with label,
    deviation, error_norm, color and optional linestyle. Units belong in label.
    The caller explicitly selects plotted models; excluded models must be
    disclosed in the companion analysis, particularly unstable truncations.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(3, len(records), figsize=(6.5 * len(records), 8.1),
                             squeeze=False, sharex='col')
    for col, record in enumerate(records):
        t = record['time']; truth = np.asarray(record['reference_deviation'])
        axes[0, col].set_title(record['label'], pad=12)
        for row, component in enumerate(coordinates):
            ax = axes[row, col]
            ax.plot(t, truth[:, component], color=RED, lw=2.8, label='Independent nonlinear ODE')
            for model in record['models']:
                ax.plot(t, np.asarray(model['deviation'])[:, component],
                        color=model['color'], ls=model.get('linestyle', '--'),
                        lw=1.3, label=model['label'])
            ax.set_ylabel(f'δx{component + 1}')
            ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useMathText=True)
            ax.grid(alpha=.2)
        for model in record['models']:
            error = np.asarray(model['error_norm'], dtype=float).copy()
            error[0] = np.nan  # Initialization is not an out-of-sample forecast.
            axes[2, col].plot(t, np.where(error > 0, error, np.nan), color=model['color'],
                              ls=model.get('linestyle', '--'), label=model['label'])
        axes[2, col].set(yscale='log', xlabel='Forcing periods since perturbation',
                         ylabel='Deviation prediction error')
        axes[2, col].grid(alpha=.2)
        axes[0, col].legend(loc='upper right', fontsize=8)
    fig.suptitle('Character-expanded transverse dynamics: free rollout after one initialization', fontsize=14)
    fig.tight_layout()
    return fig


def truncation_and_amplitude_figure(records):
    """Show supplied spectral radii and error-versus-perturbation summaries."""
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8))
    colors = [BLUE, GOLD, RED, '#009e73']
    for record, color in zip(records, colors * (1 + len(records) // len(colors))):
        axes[0].plot(record['quotients'], record['spectral_radii'], 'o-',
                     color=color, label=record['label'])
        axes[1].plot(record['perturbations'], record['rms_errors'], 'o-',
                     color=color, label=record['label'])
    axes[0].axhline(1, color='#555555', ls='--', lw=1, label='Stability boundary')
    axes[0].set(xscale='log', yscale='log', xlabel='Retained quotient period M',
                ylabel='Spectral radius over the original period N',
                title='Averaging the cocycle can create false instability')
    axes[1].set(xscale='log', yscale='log', xlabel='Initial angle displacement',
                ylabel='RMS deviation prediction error',
                title='Full cocycle: limits of the linear approximation')
    for ax in axes:
        ax.grid(alpha=.2, which='both'); ax.legend(fontsize=9)
    fig.tight_layout()
    return fig


def cyclic_partition_figure(record, *, depth=3, coordinates=(0, 1),
                            axis_labels=('x₁', 'x₂')):
    """Draw finite residue cells and native diameter/separation diagnostics.

    Requires states and geometry from cycle_partition_geometry. The terminal
    singleton level is omitted from the logarithmic curves and called out.
    This presents a finite test, not an invariant-neighborhood certificate.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap
    states = np.asarray(record['states']); geometry = record['geometry']
    levels = geometry['levels']
    if not isinstance(depth, int) or not 0 <= depth < len(levels):
        raise ValueError('Depth must be an available integer partition level')
    cells = int(levels[depth]['cells'])
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.2))
    base = ['#0072b2', '#d55e00', '#009e73', '#cc79a7', '#e69f00',
            '#56b4e9', '#343434', '#ad532d']
    cmap = ListedColormap(base[:cells]) if cells <= 8 else plt.get_cmap('turbo', cells)
    norms = BoundaryNorm(np.arange(cells + 1) - .5, cells)
    scatter = axes[0].scatter(states[:, coordinates[0]], states[:, coordinates[1]],
                             c=np.arange(len(states)) % cells, cmap=cmap, norm=norms,
                             s=48, edgecolors='white', linewidths=.4)
    axes[0].set(xlabel=axis_labels[0], ylabel=axis_labels[1],
                title=f"{record['label']}\nResidues modulo {cells}; depth {depth}")
    fig.colorbar(scatter, ax=axes[0], label='Time-index residue',
                 ticks=range(cells) if cells <= 8 else None, fraction=.045)
    nonterminal = [level for level in levels if not level['singleton_terminal_level']]
    axes[1].semilogy([v['depth'] for v in nonterminal],
                     [v['maximum_diameter'] for v in nonterminal], 'o-', color=BLUE,
                     label='Largest sampled cell diameter')
    separated = [v for v in nonterminal if v['minimum_sample_separation'] is not None]
    axes[1].semilogy([v['depth'] for v in separated],
                     [v['minimum_sample_separation'] for v in separated], 's--', color=GOLD,
                     label='Smallest sampled intercell distance')
    axes[1].set(xlabel='Binary depth n', ylabel='Euclidean distance',
                title='Shrinking finite cells — not an infinite-limit proof',
                xticks=[v['depth'] for v in nonterminal])
    axes[1].text(.03, .04, f"Depth {int(levels[-1]['depth'])}: singleton diameter = 0 automatically.\n"
                 'That terminal value is excluded.', transform=axes[1].transAxes,
                 fontsize=9, bbox=dict(facecolor='white', alpha=.9, edgecolor='#dddddd'))
    axes[1].legend(fontsize=9, loc='upper right')
    for ax in axes:
        ax.grid(alpha=.2)
    fig.tight_layout()
    return fig
