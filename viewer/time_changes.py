"""Read-only figures for native clock changes and supplied spectral references."""
from __future__ import annotations

import math


def clock_character_figure(evolution, linear_reference, *, character=0,
                           title="Characters under a change of speed"):
    """Compare a native clock query's bare/corrected character with a reference.

    linear_reference contains [real, imaginary] values on the same time mesh,
    normally from the native Fourier evolution API. This function only renders
    supplied data and discrepancies; it does not solve, fit, or invert a clock.
    """
    import matplotlib.pyplot as plt

    data = evolution.get("result", evolution)
    rows = data["samples"]
    if not isinstance(character, int) or character < 0 or not rows:
        raise ValueError("A nonnegative character index and nonempty samples are required")
    if len(linear_reference) != len(rows):
        raise ValueError("Linear reference must use the clock sample mesh")
    time = [float(row["time"]) for row in rows]
    bare = [complex(*row["bare_characters"][character]) for row in rows]
    corrected = [complex(*row["corrected_characters"][character]) for row in rows]
    reference = [complex(*value) for value in linear_reference]
    if any(b < a for a, b in zip(time, time[1:])):
        raise ValueError("Figure time values must be nondecreasing")
    if not all(math.isfinite(t) for t in time) or not all(
        math.isfinite(z.real) and math.isfinite(z.imag)
        for values in (bare, corrected, reference) for z in values
    ):
        raise ValueError("Figure data must be finite")
    fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True, constrained_layout=True)
    fig.suptitle(title, fontsize=16)
    # Contrasting colors and draw order keep an exact overlay visible.
    axes[0].plot(time, [z.real for z in bare], color="#b07800", lw=1.8,
                 label="Original character")
    axes[0].plot(time, [z.real for z in corrected], color="#d3382b", lw=2.8,
                 label="Corrected eigenfunction")
    axes[0].plot(time, [z.real for z in reference], color="#006fb0", lw=1.4,
                 ls="--", label="Constant-frequency evolution")
    axes[0].set_ylabel("Real part")
    axes[0].legend(loc="lower center", bbox_to_anchor=(.5, 1.01), ncol=3, fontsize=10)
    axes[1].plot(time, [abs(a-b) for a, b in zip(corrected, reference)],
                 color="#006fb0", lw=1.5)
    axes[1].set_ylabel("Complex prediction error")
    axes[1].set_xlabel("Physical time")
    axes[1].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    for axis in axes:
        axis.grid(alpha=.22)
        axis.set_xlim(time[0], time[-1])
    return fig
