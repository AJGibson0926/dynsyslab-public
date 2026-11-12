"""Viewer comparisons of saved phase and time traces with per-model native clocks."""
import numpy as np
from .history import matrix
from .sampling import clock
from .viewer import dashboard

def comparison_pages(cases,models,*,reference='reference',baseline=None,axes=(0,1),labels=('Value','Rate'),
                     time_label='Time',style=None,expanded=True):
    """cases[case][name]=(time,values); omit failed/empty curves and label partial coverage.
    Return a Viewer page mapping for each case plus model-scale panels, in the supplied order.
    """
    colors={'reference':('#D5252E','solid',2.3),'baseline':('#D39600','solid',4.1),'fit':('#0057BD','dash',1.7)}
    colors.update(style or {})
    pages={};models=list(models)
    if not cases or not models or len(set(models))!=len(models):raise ValueError('Expected unique models and named cases')
    for case,traces in cases.items():
        t=clock(traces[reference][0]);panels=[];zoom=[]
        for name in models:
            lines=[];fit_lines=[]
            for key,role in ([(baseline,'baseline')] if baseline is not None else [])+[(reference,'reference'),(name,'fit')]:
                ts,values=traces[key]
                if not len(ts):
                    if role=='fit':continue
                    raise ValueError('Empty reference or baseline curve')
                tt=clock(ts,1);values=matrix(values)
                if len(tt)!=len(values):raise ValueError('Unaligned comparison curve')
                label=key if tt[-1]>=t[-1] else f'{key} (to {tt[-1]:.5g})'
                x=tt if axes[0] is None else values[:,axes[0]]
                trace=(label,x,values[:,axes[1]],*colors[role]);lines.append(trace)
                if role=='fit':fit_lines.append(trace)
            panel=dict(title=f'{case} · {name}',kind='time' if axes[0] is None else 'phase',
                       x_label=time_label if axes[0] is None else labels[0],y_label=labels[1])
            panels.append(dict(panel,curves=lines))
            if fit_lines:zoom.append(dict(panel,curves=fit_lines))
        pages[case]=panels
        if expanded and zoom:pages[f'{case} · model scale']=zoom
    return pages
