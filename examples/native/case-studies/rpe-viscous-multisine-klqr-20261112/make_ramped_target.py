"""Write the ramped multi-sine RPE target as a native `trajectory` control target.

r_ref(tau) = 1 + s(tau) * (r_ms(tau) - 1), v_ref = d r_ref / d tau (exact),
with s the quintic smoothstep 0 -> 1 over [0, ramp]. r_ref is a convex blend of
the equilibrium r = 1 and r_ms >= 0.70, so it never drops below 0.70.

Usage: python make_ramped_target.py OUTPUT.npz [--ramp 10] [--stop 31] [--dt 0.001]
"""

from __future__ import annotations

import argparse
import json

import numpy as np

OFFSET = 1.5
AMPLITUDES = (0.40, 0.25, 0.15)
OMEGAS = (1.3, 2.9, 5.3)
PHASES = (0.0, 0.9, 2.1)


def multisine(t: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    r = OFFSET + sum(a * np.sin(w * t + p) for a, w, p in zip(AMPLITUDES, OMEGAS, PHASES))
    v = sum(a * w * np.cos(w * t + p) for a, w, p in zip(AMPLITUDES, OMEGAS, PHASES))
    return r, v


def smoothstep(t: np.ndarray, ramp: float) -> tuple[np.ndarray, np.ndarray]:
    x = np.clip(t / ramp, 0.0, 1.0)
    s = x**3 * (10 - 15 * x + 6 * x**2)
    ds = np.where((t > 0) & (t < ramp), 30 * x**2 * (1 - x) ** 2 / ramp, 0.0)
    return s, ds


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output")
    parser.add_argument("--ramp", type=float, default=10.0)
    parser.add_argument("--stop", type=float, default=31.0)
    parser.add_argument("--dt", type=float, default=0.001)
    args = parser.parse_args()
    t = np.round(np.arange(0.0, args.stop + 0.5 * args.dt, args.dt), 12)
    r_ms, v_ms = multisine(t)
    s, ds = smoothstep(t, args.ramp)
    values = np.column_stack([1.0 + s * (r_ms - 1.0), ds * (r_ms - 1.0) + s * v_ms])
    metadata = {
        "format": "dynsyslab.control-target-table/v1",
        "coordinates": "dimensionless (r = R/R0, v = dr/dtau, tau = t/T0)",
        "formula": "r = 1 + s(tau)(r_ms - 1); s = quintic smoothstep over [0, ramp]",
        "multisine": {"offset": OFFSET, "amplitudes": AMPLITUDES,
                      "angular_frequencies": OMEGAS, "phases_rad": PHASES},
        "ramp": args.ramp, "r_min": float(values[:, 0].min()), "r_max": float(values[:, 0].max()),
    }
    np.savez(args.output, time=t, control_target_values=values,
             metadata_json=np.array(json.dumps(metadata)))
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
