"""Export chronological, error-bounded vertices of any two coordinates plus their source clock."""
import argparse
import json
from pathlib import Path
import tempfile
import numpy as np
from resolution import DisplayVertices, read_trace


def export(source, output, library, axes, time_key='time', state_key='values', tolerance=2e-5):
    t,x=read_trace(source,time_key,state_key)
    if t.ndim!=1 or len(t)!=len(x) or not np.isfinite(t).all() or np.any(np.diff(t)<=0):
        raise ValueError('Expected a finite, increasing sample clock')
    ix=DisplayVertices(library)(x[:,axes[0]],x[:,axes[1]],tolerance)
    if output.exists():
        raise FileExistsError(output)
    meta=dict(source=str(source.resolve()),axes=list(axes),sample_count=len(t),display_count=len(ix),
              tolerance=tolerance,time_key=time_key,state_key=state_key)
    output.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(dir=output.parent) as tmp:
        ready=Path(tmp)/'trace.npz'
        np.savez_compressed(ready,time=t[ix],state=x[ix][:,axes],source_rows=ix,metadata_json=json.dumps(meta))
        ready.rename(output)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('source','out','library'):p.add_argument('--'+key,required=True,type=Path)
    p.add_argument('--axes',type=int,nargs=2,default=[0,1])
    p.add_argument('--time-key',default='time')
    p.add_argument('--state-key',default='values')
    p.add_argument('--tolerance',type=float,default=2e-5)
    a=p.parse_args()
    export(a.source,a.out,a.library,a.axes,a.time_key,a.state_key,a.tolerance)
