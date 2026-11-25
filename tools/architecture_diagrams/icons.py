"""Conceptual icons for architecture-diagram blocks.

Every icon draws into a 100 x 100 box and is returned as SVG markup. Add an icon
by writing one function and registering it in ICONS.
"""
from __future__ import annotations

import math

BLUE, RED, GOLD = "#1F5FBF", "#D8242A", "#D9A21B"
INK, MUTED, PAPER = "#26323F", "#7A8794", "#FFFFFF"


def _f(value: float) -> str:
    return f"{value:.2f}".rstrip("0").rstrip(".")


def _polyline(points, colour, width, extra="") -> str:
    path = " ".join(f"{_f(x)},{_f(y)}" for x, y in points)
    return (f'<polyline points="{path}" fill="none" stroke="{colour}" stroke-width="{width}" '
            f'stroke-linecap="round" stroke-linejoin="round" {extra}/>')


def monitor() -> str:
    """Small computer monitor with a plot on screen (solver, post-processor)."""
    screen = _polyline([(18 + 3 * i, 46 - 14 * math.sin(i / 2.2) * math.exp(-i / 14)) for i in range(22)], BLUE, 2.6)
    return (f'<rect x="8" y="12" width="84" height="60" rx="5" fill="{INK}"/>'
            f'<rect x="13" y="17" width="74" height="50" rx="2" fill="#EAF1FB"/>{screen}'
            f'<rect x="43" y="72" width="14" height="10" fill="{MUTED}"/>'
            f'<rect x="28" y="82" width="44" height="6" rx="3" fill="{INK}"/>')


def matrix() -> str:
    """Hankel matrix: bracketed grid whose anti-diagonals share a shade (H-DMD)."""
    shades = ["#DCE7F7", "#B9CFEF", "#8FB2E4", "#5F8FD3", "#3A6FC4", BLUE]
    cells = []
    for row in range(6):
        for col in range(6):
            colour = shades[(row + col) % len(shades)]
            cells.append(f'<rect x="{20 + 10 * col}" y="{20 + 10 * row}" width="9" height="9" fill="{colour}"/>')
    bracket = (f'<path d="M18 16 H12 V84 H18 M82 16 H88 V84 H82" fill="none" stroke="{INK}" '
               f'stroke-width="3" stroke-linecap="round"/>')
    return "".join(cells) + bracket


def library() -> str:
    """Library shelf of candidate functions (SINDy)."""
    books = [(14, 30, BLUE), (24, 22, RED), (33, 34, GOLD), (43, 26, "#5F8FD3"), (52, 30, INK),
             (61, 20, RED), (70, 32, BLUE), (79, 24, GOLD)]
    shelf = [f'<rect x="8" y="{y}" width="84" height="4" fill="{INK}"/>' for y in (48, 88)]
    out = [f'<rect x="8" y="10" width="4" height="82" fill="{INK}"/>',
           f'<rect x="88" y="10" width="4" height="82" fill="{INK}"/>'] + shelf
    for x, height, colour in books:
        out.append(f'<rect x="{x}" y="{48 - height}" width="8" height="{height}" rx="1" fill="{colour}"/>')
    for i, (x, height, colour) in enumerate(reversed(books)):
        out.append(f'<rect x="{x}" y="{88 - height + 4 * (i % 2)}" width="8" height="{height - 4 * (i % 2)}" '
                   f'rx="1" fill="{colour}" opacity="0.85"/>')
    return "".join(out)


def spectrum() -> str:
    """Koopman spectrum: dashed unit circle in the complex plane, red eigenvalues."""
    # Conjugate pairs on and inside the unit circle at many angles (a typical H-DMD spectrum).
    polar = [(0.99, 0.0), (0.98, 0.35), (0.97, 0.7), (0.99, 1.05), (0.96, 1.45), (0.98, 1.85), (0.95, 2.3),
             (0.97, 2.75), (0.99, 3.1), (0.8, 0.55), (0.72, 1.25), (0.84, 2.05), (0.62, 2.6), (0.55, 0.15),
             (0.4, 1.7), (0.3, 2.9)]
    dots = []
    for radius, angle in polar:
        for sign in ((1,) if angle in (0.0, 3.1) else (1, -1)):
            dots.append(f'<circle cx="{_f(50 + 36 * radius * math.cos(angle))}" '
                        f'cy="{_f(50 - sign * 36 * radius * math.sin(angle))}" r="2.6" fill="{RED}"/>')
    return (f'<line x1="6" y1="50" x2="94" y2="50" stroke="{MUTED}" stroke-width="1.4"/>'
            f'<line x1="50" y1="6" x2="50" y2="94" stroke="{MUTED}" stroke-width="1.4"/>'
            f'<circle cx="50" cy="50" r="36" fill="none" stroke="{INK}" stroke-width="2" stroke-dasharray="5 4"/>'
            + "".join(dots))


def _diverging(value: float) -> str:
    value = max(-1.0, min(1.0, value))
    low, mid, high = (31, 95, 191), (247, 247, 247), (216, 36, 42)
    a, b = (mid, high) if value >= 0 else (mid, low)
    t = abs(value)
    rgb = tuple(round(a[i] + (b[i] - a[i]) * t) for i in range(3))
    return "#%02X%02X%02X" % rgb


def heatmap() -> str:
    """Eigenfunction on a phase plane: Re psi as a spiralling diverging heat map."""
    n, cells = 16, []
    for row in range(n):
        for col in range(n):
            x, y = (col + 0.5) / n * 2 - 1, 1 - (row + 0.5) / n * 2
            r = math.hypot(x, y)
            value = math.cos(math.atan2(y, x) + 2.4 * r) * min(1.0, 1.4 * r)
            cells.append(f'<rect x="{_f(10 + 80 * col / n)}" y="{_f(10 + 80 * row / n)}" '
                         f'width="{_f(80 / n + 0.3)}" height="{_f(80 / n + 0.3)}" fill="{_diverging(value)}"/>')
    return "".join(cells) + f'<rect x="10" y="10" width="80" height="80" fill="none" stroke="{INK}" stroke-width="1.6"/>'


def time_series() -> str:
    """Viewer-style trajectory plot with several initial conditions overlaid."""
    curves = []
    for amplitude, phase, decay, colour in ((26, 0.0, 22, BLUE), (18, 1.4, 30, RED), (12, 2.6, 40, GOLD)):
        curves.append(_polyline([(14 + 1.5 * i, 50 - amplitude * math.cos(i / 4.0 + phase) * math.exp(-i / decay))
                                 for i in range(52)], colour, 2.4))
    return (f'<rect x="6" y="8" width="88" height="84" fill="{PAPER}"/>'
            f'<polyline points="12,12 12,88 92,88" fill="none" stroke="{INK}" stroke-width="2"/>' + "".join(curves))


def krom() -> str:
    """K-ROM: diagonal Lambda, input column B_psi and observable rows C."""
    diag = "".join(f'<rect x="{12 + 11 * i}" y="{12 + 11 * i}" width="10" height="10" fill="{RED}"/>' for i in range(5))
    lam = f'<rect x="10" y="10" width="58" height="58" fill="#FBE9EA" stroke="{INK}" stroke-width="1.6"/>{diag}'
    b = "".join(f'<rect x="78" y="{12 + 11 * i}" width="10" height="10" fill="{BLUE}"/>' for i in range(5))
    b = f'<rect x="76" y="10" width="14" height="58" fill="#E6EEF9" stroke="{INK}" stroke-width="1.6"/>{b}'
    c = "".join(f'<rect x="{12 + 11 * i}" y="{78 + 9 * row}" width="10" height="7" fill="{GOLD}" '
                f'opacity="{0.55 + 0.45 * ((i + row) % 2)}"/>' for i in range(5) for row in range(2))
    c = f'<rect x="10" y="76" width="58" height="20" fill="#FBF3DF" stroke="{INK}" stroke-width="1.6"/>{c}'
    label = (f'<text x="39" y="9" font-size="9" font-style="italic" text-anchor="middle" fill="{INK}">Λ</text>'
             f'<text x="83" y="9" font-size="9" font-style="italic" text-anchor="middle" fill="{INK}">B</text>'
             f'<text x="96" y="90" font-size="9" font-style="italic" text-anchor="middle" fill="{INK}">C</text>')
    return lam + b + c + label


def controller() -> str:
    """Feedback loop: reference, summing junction, gain K, plant, measured return."""
    arrow = 'marker-end="url(#icon-arrow)"'
    return (f'<circle cx="22" cy="40" r="8" fill="{PAPER}" stroke="{INK}" stroke-width="2"/>'
            f'<path d="M17 40 H27 M22 35 V45" stroke="{INK}" stroke-width="1.6"/>'
            f'<rect x="38" y="28" width="20" height="24" rx="3" fill="{BLUE}"/>'
            f'<text x="48" y="45" font-size="14" font-weight="700" text-anchor="middle" fill="{PAPER}">K</text>'
            f'<rect x="66" y="28" width="26" height="24" rx="3" fill="{INK}"/>'
            f'<path d="M70 44 Q75 32 79 40 T88 36" fill="none" stroke="{GOLD}" stroke-width="2"/>'
            f'<path d="M4 40 H12" stroke="{INK}" stroke-width="2" {arrow}/>'
            f'<path d="M30 40 H36" stroke="{INK}" stroke-width="2" {arrow}/>'
            f'<path d="M58 40 H64" stroke="{INK}" stroke-width="2" {arrow}/>'
            f'<path d="M92 40 H96 V74 H22 V50" fill="none" stroke="{RED}" stroke-width="2" {arrow}/>')


def estimator() -> str:
    """Magnifier recovering a smooth hidden state from scattered measurements."""
    dots = "".join(f'<circle cx="{_f(10 + 6 * i)}" cy="{_f(52 - 22 * math.sin(i / 2.0) + (5 if i % 3 else -4))}" '
                   f'r="2.2" fill="{MUTED}"/>' for i in range(14))
    curve = _polyline([(10 + 1.5 * i, 52 - 22 * math.sin(i / 8.0)) for i in range(56)], BLUE, 2.6)
    return (dots + curve + f'<circle cx="62" cy="48" r="20" fill="#FFFFFF" fill-opacity="0.35" '
            f'stroke="{INK}" stroke-width="3"/>'
            f'<path d="M76 62 L92 80" stroke="{INK}" stroke-width="6" stroke-linecap="round"/>')


def archive() -> str:
    """Archive record: a stack of .npz files (arrays + metadata)."""
    out = []
    for i, colour in enumerate(("#B9CFEF", "#8FB2E4", BLUE)):
        x, y = 22 + 6 * i, 14 + 6 * i
        out.append(f'<path d="M{x} {y} h38 l12 12 v52 h-50 z" fill="{colour}" stroke="{INK}" stroke-width="2"/>')
    out.append(f'<path d="M72 26 h-12 v-12" fill="none" stroke="{INK}" stroke-width="2"/>')
    out.append(f'<text x="59" y="66" font-size="13" font-weight="700" text-anchor="middle" fill="{PAPER}">.npz</text>')
    return "".join(out)


def svd() -> str:
    """Truncated SVD: U (tall), Sigma (diagonal), V^T (wide)."""
    diag = "".join(f'<rect x="{41 + 6 * i}" y="{31 + 6 * i}" width="5" height="5" fill="{RED}"/>' for i in range(4))
    return (f'<rect x="6" y="20" width="24" height="60" fill="#DCE7F7" stroke="{INK}" stroke-width="2"/>'
            f'<rect x="38" y="28" width="28" height="28" fill="#FBE9EA" stroke="{INK}" stroke-width="2"/>{diag}'
            f'<rect x="72" y="28" width="22" height="28" fill="#FBF3DF" stroke="{INK}" stroke-width="2"/>'
            f'<text x="18" y="94" font-size="12" font-style="italic" text-anchor="middle" fill="{INK}">U</text>'
            f'<text x="52" y="72" font-size="12" font-style="italic" text-anchor="middle" fill="{INK}">Σ</text>'
            f'<text x="83" y="72" font-size="12" font-style="italic" text-anchor="middle" fill="{INK}">Vᵀ</text>')


ICONS = {"archive": archive, "svd": svd, "monitor": monitor, "matrix": matrix, "library": library, "spectrum": spectrum,
         "heatmap": heatmap, "time_series": time_series, "krom": krom, "controller": controller,
         "estimator": estimator}


def icon_markup(name: str, x: float, y: float, size: float) -> str:
    scale = size / 100.0
    return f'<g transform="translate({_f(x)},{_f(y)}) scale({_f(scale)})">{ICONS[name]()}</g>'


ICON_DEFS = (f'<marker id="icon-arrow" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="5" markerHeight="5" '
             f'orient="auto-start-reverse"><path d="M0 0 L10 5 L0 10 z" fill="{INK}"/></marker>')
