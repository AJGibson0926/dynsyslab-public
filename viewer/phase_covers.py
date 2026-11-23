"""Presentation-only phase-cover explainer for precomputed planar orbit samples."""
import json
import math
import re
from pathlib import Path


def _validate_phase_geometries(geometries):
    if not geometries:raise ValueError('Supply phase-cover geometries')
    for key,data in geometries.items():
        if data.get('schema')!='dynsyslab.phase-cover-geometry/v1' or str(data['prime'])!=str(key):raise ValueError('Invalid phase-cover schema')
        if data['prime'] not in (2,3) or len(data['orbits'])!=3:raise ValueError('This explainer displays three levels for prime 2 or 3')
        if not isinstance(data['samples'],int) or data['samples']<3:raise ValueError('Supply at least three samples per orbit')
        for orbit in data['orbits']:
            if orbit['state_dimension']!=2 or len(orbit['values'])!=2*data['samples']:raise ValueError('Planar sampled geometry required')


def _presentation(geometries, axis_labels, bounds):
    """Resolve display extents only; no physical coordinate map is inferred."""
    if len(axis_labels)!=2 or any(not isinstance(s,str) or not s for s in axis_labels):
        raise ValueError('Supply two nonempty axis labels')
    if bounds is None:
        values=[orbit['values'] for data in geometries.values() for orbit in data['orbits']]
        bounds=[]
        for j in range(2):
            lo=min(min(v[j::2]) for v in values);hi=max(max(v[j::2]) for v in values)
            pad=.06*(hi-lo) if hi>lo else max(1,abs(lo))*.06
            bounds.extend((lo-pad,hi+pad))
    if len(bounds)!=4 or any(not math.isfinite(v) for v in bounds) or bounds[0]>=bounds[1] or bounds[2]>=bounds[3]:
        raise ValueError('Bounds must be finite [xmin,xmax,ymin,ymax] with positive spans')
    return dict(axisLabels=list(axis_labels),bounds=list(bounds))


def phase_cover_fragment(geometries, *, root_id='phase-cover-view',
                        axis_labels=('x₁','x₂'), bounds=None):
    """Render supplied planar covers with time, phase, zoom and drag controls.

    Uses sample_phase_cover's output directly; optional historical branch tables
    are ignored. Display bounds default to the combined supplied geometry.
    """
    if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_-]*',root_id):raise ValueError('Invalid root id')
    _validate_phase_geometries(geometries)
    presentation=json.dumps(_presentation(geometries,axis_labels,bounds),ensure_ascii=False).replace('<','\\u003c')
    payload=json.dumps(geometries,allow_nan=False,separators=(',',':')).replace('<','\\u003c')
    template=Path(__file__).with_name('phase_cover_fragment.html').read_text(encoding='utf-8')
    clock=Path(__file__).with_name('phase_cover_clock.js').read_text(encoding='utf-8')
    return (template.replace('__ROOT_ID__',root_id).replace('__GEOMETRY__',payload)
            .replace('__CLOCK_HELPER__',clock).replace('__PRESENTATION__',presentation))


def phase_cover_figure(geometries, *, prime=2, base_phase_turns=1/6, digits=(0, 1),
                       axis_labels=('x₁','x₂'), bounds=None):
    """Static notebook/export counterpart using the same precomputed geometry.

    Markers interpolate sampled curves for display only. Colors show mapped
    base phase, not geometric polar angle. The arithmetic is finite labeling;
    no ODE or inverse frequency calculation occurs here.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    _validate_phase_geometries(geometries)
    display=_presentation(geometries,axis_labels,bounds)
    data = geometries[str(prime)]
    if len(digits) != 2 or any(not isinstance(a, int) or not 0 <= a < prime for a in digits):
        raise ValueError('Supply two base-p digits')
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.7), sharex=True, sharey=True)
    palette = np.array([[.84, .15, .16], [0, .55, .35], [0, .45, .70]])
    for j, (ax, orbit) in enumerate(zip(axes, data['orbits'])):
        x = np.asarray(orbit['values']).reshape(-1, 2)
        n = len(x); count = prime ** j
        segments = np.stack([x, np.roll(x, -1, axis=0)], axis=1)
        u = np.mod(count * (np.arange(n) + .5) / n, 1) * 3
        k = np.floor(u).astype(int); blend = (u - k)[:, None]
        colors = (1 - blend) * palette[k] + blend * palette[(k + 1) % 3]
        ax.add_collection(LineCollection(segments, colors=colors, linewidths=2.2))
        selected = 0 if j == 0 else digits[0] + (prime * digits[1] if j == 2 else 0)
        for residue in range(count):
            turn = ((base_phase_turns + residue) / count) % 1
            index = turn * n; low = int(index); frac = index - low
            point = (1 - frac) * x[low] + frac * x[(low + 1) % n]
            ax.scatter(*point, s=70, facecolor='#17232e' if residue == selected else 'white',
                       edgecolor='#17232e', linewidth=1.5, zorder=3)
            ax.annotate(str(residue), point, xytext=(6, 7), textcoords='offset points')
        ax.set(xlim=display['bounds'][:2], ylim=display['bounds'][2:], xlabel=axis_labels[0],
               title=f'ω₀/{count}: {count} compatible phase point(s)')
        ax.grid(alpha=.2)
    axes[0].set_ylabel(axis_labels[1])
    fig.suptitle(f'Chosen phase-cover tower, p={prime}; base phase={360*base_phase_turns:g}°')
    fig.tight_layout()
    return fig


def phase_cover_evolution_figure(geometries, *, prime=2, initial_phase_turns=1/6,
                                digits=(1,1), elapsed_base_periods=0,
                                axis_labels=('x₁','x₂'), bounds=None):
    """Animate the exact clock on a finite displayed prefix, with p-adic carries.

    This only relabels/interpolates supplied geometry; it does not infer a tower,
    integrate an ODE, or certify infinite dynamics. Time is in base-orbit periods.
    """
    import math
    if prime not in (2,3) or len(digits)!=2 or any(not isinstance(a,int) or not 0<=a<prime for a in digits):
        raise ValueError('Provide a supported prime and two valid digits')
    total=initial_phase_turns+elapsed_base_periods
    if not math.isfinite(total):raise ValueError('Finite clock required')
    turns=math.floor(total)
    residue=(digits[0]+prime*digits[1]+turns)%(prime**2)
    current=(residue%prime,residue//prime)
    fig=phase_cover_figure(geometries,prime=prime,base_phase_turns=total-turns,digits=current,
                          axis_labels=axis_labels,bounds=bounds)
    fig.suptitle(f'Solenoid flow: t/T₀={elapsed_base_periods:g}; residues modulo {prime}, {prime**2}: {current[0]}, {residue}')
    return fig


def solenoid_suspension_figure(*, base=2, depth=3, phase_turns=.4, address=3):
    """Cut-open suspension schematic and finite Cantor-fiber construction.

    This displays a classical adding-machine suspension, not a physical phase
    portrait or a certified attractor. Each row stands for a cylinder of
    unspecified higher digits. Vertical placement is a schematic Cantor
    embedding; it is not energy, radius, or a Euclidean metric on the system.
    """
    import matplotlib.pyplot as plt
    if type(base) is not int or not 2 <= base <= 32 or type(depth) is not int or not 1 <= depth <= 5 or base**depth > 32:
        raise ValueError('Require base>=2, depth>=1 and at most 32 displayed prefixes')
    if not math.isfinite(phase_turns) or not 0 <= phase_turns < 1 or type(address) is not int:
        raise ValueError('Supply phase in [0,1) and an integer address')
    n=base**depth; radix=2*base-1
    def interval(a,d):
        y=sum(2*((a//base**j)%base)/radix**(j+1) for j in range(d))
        return y,y+radix**(-d)
    ys=[sum(interval(a,depth))/2 for a in range(n)]
    chosen=address%n; successor=(address+1)%n
    fig,(ax,fiber)=plt.subplots(1,2,figsize=(12,6.2),width_ratios=(3,1.15),sharey=True)
    for a,y in enumerate(ys):
        ax.plot([0,1],[y,y],color='#aeb9bf',lw=1)
        ax.scatter([phase_turns],[y],color='#c52b24',s=32,zorder=4)
    ax.axvline(phase_turns,color='#c52b24',lw=1,ls=':',zorder=1)
    ax.axvline(0,color='#7c8589',lw=1,ls='--')
    ax.axvline(1,color='#7c8589',lw=1,ls='--')
    for lo,hi,a in ((phase_turns,1,chosen),(.02,.22,successor)):
        ax.annotate('',xy=(hi,ys[a]),xytext=(lo,ys[a]),
                    arrowprops={'arrowstyle':'->','color':'#0072b2','lw':2.5})
    ax.scatter([1,0],[ys[chosen],ys[successor]],color='#0072b2',s=65,zorder=5)
    labeled=list(range(n)) if n<=8 else sorted({chosen,successor})
    ax.set(xlim=(-.025,1.025),ylim=(-.09,1.08),xlabel='Base phase / 2π — physical time runs to the right',
           ylabel=f'Address modulo {n} (schematic transverse placement)',
           yticks=[ys[a] for a in labeled],yticklabels=[str(a) for a in labeled],
           title='Flow segments; the two edges are glued')
    ax.text(phase_turns,1.015,'One fiber',color='#c52b24',ha='center')
    for d in range(1,depth+2):
        for a in range(base**d):
            lo,hi=interval(a,d)
            fiber.plot([d,d],[lo,hi],color='#c52b24',lw=7,solid_capstyle='butt')
    fiber.set(xlim=(.5,depth+1.65),xticks=range(1,depth+2),
              xlabel='Refinement depth',title='Fiber: Cantor refinement')
    fiber.tick_params(axis='y',left=False,labelleft=False)
    fiber.text(depth+1.45,.5,'⋮',ha='center',va='center')
    for axis in (ax,fiber):
        axis.spines[['top','right']].set_visible(False)
    fig.suptitle(f'{base}-adic solenoid, cut open — finite-prefix schematic',fontsize=16)
    fig.text(.5,.055,f'Seam: right edge at address a={address} continues at left edge a={address+1}. '
             'Higher digits remain unspecified.',ha='center',fontsize=11)
    fig.text(.5,.018,'The red dotted guide is not a path: the fiber contains Cantor points, with gaps at every depth.',
             ha='center',fontsize=11)
    fig.tight_layout(rect=(0,.10,1,.94))
    return fig
