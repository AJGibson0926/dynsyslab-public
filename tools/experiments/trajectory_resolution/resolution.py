"""Bound display error for 2-D paths and audit time-sampled traces from any system."""
import argparse
import ctypes
import json
from pathlib import Path

import numpy as np


class DisplayVertices:
    def __init__(self, library):
        self.library = ctypes.CDLL(str(Path(library).resolve()))
        self.select = self.library.trajectory_vertices
        pointer = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='ALIGNED')
        self.select.argtypes = [pointer, ctypes.c_ssize_t, pointer, ctypes.c_ssize_t,
                                ctypes.c_size_t, ctypes.c_double,
                                np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags='C_CONTIGUOUS')]
        self.select.restype = ctypes.c_size_t

    def __call__(self, x, y, tolerance=2e-5):
        x, y = (np.require(a, dtype='float64', requirements='A') for a in (x,y))
        if x.ndim!=1 or y.shape!=x.shape or not len(x):
            raise ValueError('Expected two equally long, nonempty vectors')
        indices = np.empty(len(x),dtype=np.int64)
        count = self.select(x, x.strides[0]//8, y, y.strides[0]//8, len(x), tolerance, indices)
        if not count:
            raise ValueError('Expected finite coordinates and a fractional tolerance between zero and one')
        return indices[:count].copy()


def read_trace(path, time_key='time', state_key='values'):
    with np.load(path, allow_pickle=False) as data:
        t, x = data[time_key], data[state_key]
    return t, x[:,None] if x.ndim==1 else x


def compare_clocks(reference_time, reference, sample_time, samples, window=None):
    """Compare a sampled polyline to reference values at common times without extrapolation."""
    rt, st = (np.asarray(t, dtype=float) for t in (reference_time, sample_time))
    ref, sampled = (np.asarray(a, dtype=float) for a in (reference,samples))
    if any(t.ndim!=1 or len(t)<2 or not np.isfinite(t).all() or np.any(np.diff(t)<=0) for t in (rt,st)):
        raise ValueError('Each clock must be finite and strictly increasing')
    if (ref.ndim!=2 or sampled.ndim!=2 or ref.shape[0]!=len(rt) or sampled.shape!=(len(st),ref.shape[1])
            or not np.isfinite(ref).all() or not np.isfinite(sampled).all()):
        raise ValueError('Expected finite sample-by-coordinate arrays on both clocks')
    lo, hi = max(rt[0],st[0]), min(rt[-1],st[-1])
    if window is not None:
        if len(window)!=2 or not np.isfinite(window).all() or window[1]<window[0]:
            raise ValueError('Invalid time window')
        lo,hi=max(lo,window[0]),min(hi,window[1])
    selection = slice(rt.searchsorted(lo),rt.searchsorted(hi,side='right'))
    t, x = rt[selection],ref[selection]
    if not len(t):
        raise ValueError('No reference samples in the common window')
    error=[np.interp(t,st,component)-values for component,values in zip(sampled.T,x.T)]
    return dict(compared_time=[float(t[0]),float(t[-1])], reference_count=len(t),
                sample_count=int(st.searchsorted(hi,side='right')-st.searchsorted(lo)),
                max_abs=[float(abs(e).max()) for e in error],
                rms=[float(np.sqrt(np.mean(e*e))) for e in error],
                worst_time=[float(t[abs(e).argmax()]) for e in error])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('reference',type=Path)
    parser.add_argument('sampled',type=Path)
    for which in ('reference','sample'):
        parser.add_argument('--'+which+'-time', default='time')
        parser.add_argument('--'+which+'-values', default='values')
    parser.add_argument('--window',nargs=2,type=float)
    args=parser.parse_args()
    result=compare_clocks(*read_trace(args.reference,args.reference_time,args.reference_values),
                          *read_trace(args.sampled,args.sample_time,args.sample_values),args.window)
    print(json.dumps(result,indent=2))
