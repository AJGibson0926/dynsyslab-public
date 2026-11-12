"""Match the source trajectory, the receiver clock and any finer replay windows."""
from dataclasses import dataclass
import json
from pathlib import Path
import sys

import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'trajectory_resolution'))
from resolution import DisplayVertices, read_trace
from study import TERMS, measurements


@dataclass
class PhaseWindow:
    time: np.ndarray
    state: np.ndarray
    query_time: np.ndarray
    prediction: np.ndarray
    leading: np.ndarray
    total: np.ndarray


def interval(clock,start,stop):
    return slice(clock.searchsorted(start),clock.searchsorted(stop,side='right'))


def phase_windows(folder,held,prediction,bounds,nshow,replays):
    native_t,source=read_trace(held['meta']['source'])
    t=held['time']
    total,_=measurements(held,range(len(TERMS)))
    total/=held['ranges_m']*1000.
    mono=held['monopole_pressure']/1000.
    windows=[]
    for i,(lo,hi) in enumerate([(max(t[0],bounds[0]),bounds[nshow]),
                                (max(t[0],bounds[-nshow-1]),bounds[-1]),(t[0],t[-1])]):
        idx=interval(native_t,lo,hi)
        nt,xx=native_t[idx],source[idx]
        if i<2 and replays:
            with np.load(replays[i]) as z:
                meta=json.loads(z['metadata_json'].item())
                if (meta['schema']!='dynsyslab.receiver-replay/v1' or meta['case']!='all'
                    or Path(meta['study']).resolve()!=folder.resolve() or meta['source']!=held['meta']['source']):
                    raise ValueError('Replay must match this study, source and correction case')
                if not np.array_equal(nt,z['time']) or not np.array_equal(xx,z['state']):
                    raise ValueError('Replay must cover the requested time window at native resolution')
                windows.append(PhaseWindow(nt,xx,nt,z['prediction'],z['monopole']/1000.,z['pressure']/1000.))
        else:
            idx=interval(t,lo,hi)
            windows.append(PhaseWindow(nt,xx,t[idx],prediction[idx],mono[idx],total[idx]))
    return windows
