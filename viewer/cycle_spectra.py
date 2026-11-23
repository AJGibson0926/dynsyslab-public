"""Presentation of precomputed finite-cycle spectra and predictions.

No integrations, decompositions, cycle detection, or spectral estimation live
here. Matplotlib is imported lazily for notebooks and exportable figures.
Records contain a label and optional eigenvalues, atoms, measures, states,
truth, and predictions. Complex eigenvalues are [real, imaginary] pairs.
"""
from __future__ import annotations

COLORS = {"resdmd": "#c52b24", "mpedmd": "#0072b2"}


def cycle_clock_fragment(periods=(2, 4, 8), *, root_id="cycle-clock",
                         step_label="Advance one step", initial_step=6):
    """Present finite quotient clocks of an adding machine, not physical orbits.

    Periods form a divisibility chain. The browser counts with exact integers;
    circular positions only visualize residues. This does not infer a cyclic
    partition, conjugacy, or infinite limit from a dynamical system's data.
    """
    import json
    import re
    from pathlib import Path
    periods = tuple(periods)
    if not 1 <= len(periods) <= 6 or any(type(n) is not int or not 2 <= n <= 32 for n in periods):
        raise ValueError("Supply one to six integer periods between 2 and 32")
    if any(b <= a or b % a for a, b in zip(periods, periods[1:])):
        raise ValueError("Periods must be a strictly increasing divisibility chain")
    if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_-]*", root_id):
        raise ValueError("Invalid root id")
    if type(initial_step) is not int or not isinstance(step_label, str) or not step_label:
        raise ValueError("Supply an integer initial step and a nonempty step label")
    payload = json.dumps(dict(periods=periods, initialStep=str(initial_step),
                              stepLabel=step_label), ensure_ascii=True).replace("<", "\\u003c")
    template = Path(__file__).with_name("cycle_clock_fragment.html").read_text(encoding="utf-8")
    return template.replace("__ROOT_ID__", root_id).replace("__CONFIG__", payload)


def spectrum_sequence(records, *, title="Finite-cycle Koopman eigenvalues"):
    """Show fitted eigenangles alongside independently supplied reference atoms."""
    import numpy as np
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(12, max(4, .8 * len(records))))
    for j, record in enumerate(records):
        angles = np.asarray(record.get("atoms", [])) / (2 * np.pi)
        ax.scatter(angles, np.full(len(angles), j), marker="|", s=110,
                   color="#737373", label="Finite-cycle reference" if j == 0 else None)
        for offset, (method, pairs) in zip([-.14, .14], record["eigenvalues"].items()):
            z = np.asarray(pairs)
            angle = np.mod(np.arctan2(z[:, 1], z[:, 0]) / (2 * np.pi), 1)
            # Put values numerically indistinguishable from one at zero.
            angle = np.where(abs(angle - 1) < 1e-8, 0, angle)
            ax.scatter(angle, np.full(len(angle), j + offset), s=15,
                       marker="x" if method == "resdmd" else "+",
                       color=COLORS.get(method, "#b57b00"), label=method if j == 0 else None)
    ax.set(yticks=range(len(records)), yticklabels=[r["label"] for r in records],
           xlabel="Eigenangle / (2π), per return-map step", xlim=(-.025, 1.025), title=title)
    ax.grid(axis="x", alpha=.2)
    ax.legend(loc="upper left", bbox_to_anchor=(0, 1.16), ncol=3)
    fig.tight_layout()
    return fig


def return_sections(records, *, coordinates=(0, 1), labels=("x₁", "x₂")):
    """Draw supplied return-map states, colored by chronological sample index."""
    import numpy as np
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, len(records), figsize=(4 * len(records), 4), squeeze=False)
    for ax, record in zip(axes[0], records):
        x = np.asarray(record["states"])
        ax.scatter(x[:, coordinates[0]], x[:, coordinates[1]], c=range(len(x)),
                   cmap="turbo", s=25, edgecolors="white", linewidths=.3)
        ax.set(xlabel=labels[0], ylabel=labels[1], title=record["label"])
        ax.grid(alpha=.2)
    fig.suptitle("Stroboscopic states — one sample per forcing period")
    fig.tight_layout()
    return fig


def prediction_comparison(record, *, samples=128):
    """Show precomputed free predictions and their pointwise differences."""
    import numpy as np
    import matplotlib.pyplot as plt
    truth = np.asarray(record["truth"])
    n = min(samples, len(truth))
    fig, axes = plt.subplots(2, 1, figsize=(12, 6), sharex=True,
                             gridspec_kw={"height_ratios": [2, 1]})
    axes[0].plot(truth[:n], color="#b57b00", lw=3, label="Held-out trajectory")
    for method, values in record["predictions"].items():
        values = np.asarray(values)
        if len(values) != len(truth):
            raise ValueError("Predictions and truth must have equal lengths")
        axes[0].plot(values[:n], color=COLORS.get(method), ls="--", lw=1.3, label=method)
        axes[1].plot((values - truth)[:n], color=COLORS.get(method), label=method)
    axes[0].set_ylabel(record.get("observable", "Observable"))
    axes[0].set_title(record["label"], pad=44)
    axes[0].legend(ncol=3, loc="lower left", bbox_to_anchor=(0, 1.01))
    axes[1].set(xlabel="Return-map steps from initialization", ylabel="Prediction − truth")
    for ax in axes:
        ax.grid(alpha=.2)
    fig.tight_layout()
    return fig


def regularized_measures(records, *, linthresh=1e-5):
    """Plot supplied signed densities; retain negative lobes and absolute mass."""
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import SymmetricalLogLocator
    fig, axes = plt.subplots(len(records), 1, figsize=(12, 2.5 * len(records)),
                             sharex=True, squeeze=False)
    for ax, record in zip(axes[:, 0], records):
        for method, measure in record["measures"].items():
            ax.plot(np.asarray(measure["grid"]) / (2 * np.pi), measure["density"],
                    color=COLORS.get(method, "#0072b2"), lw=1.1, label=method)
        for angle in record.get("atoms", []):
            t = (angle + np.pi) % (2 * np.pi) - np.pi
            ax.axvline(t / (2 * np.pi), color="#777777", alpha=.25, lw=.6)
        ax.set(title=record["label"], ylabel="Signed density")
        ax.set_yscale("symlog", linthresh=linthresh)
        ax.yaxis.set_major_locator(SymmetricalLogLocator(base=100, linthresh=linthresh))
        ax.grid(alpha=.15)
    axes[-1, 0].set_xlabel("Spectral angle / (2π), per return-map step")
    fig.suptitle("Rigged spectral regularization of mpEDMD; gray lines: cycle atoms")
    fig.tight_layout()
    return fig
