"""Run a standalone rank sweep from labeled observation records."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from miscellaneous_computes.history import fit_basis,HistoryBasis
from miscellaneous_computes.study import rank_sweep
from miscellaneous_computes.controls import run_controls

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('manifest',type=Path);ap.add_argument('output',type=Path)
    ap.add_argument('--library',type=Path,required=True,help='Compiled misc_kernels library')
    args=ap.parse_args()
    root=args.manifest.resolve().parent;cfg=json.loads(args.manifest.read_text(encoding='utf-8'))
    if cfg.get('schema')!='dynsyslab.miscellaneous-computes/study/v1':raise ValueError('Unsupported study schema')
    def load(name):
        with np.load(root/name,allow_pickle=False) as z:
            y=z['values'];y=y[:,None] if y.ndim==1 else y
            v=z['rates'];v=v[:,None] if v.ndim==1 else v
            return z['time'],y,v
    splits={role:cfg[role] for role in ('training','validation','test')}
    paths=[(root/p).resolve() for group in splits.values() for p in group.values()]
    if len(set(paths))!=len(paths):raise ValueError('A trajectory path is repeated across the fit/evaluation splits')
    by_id={k:load(p) for group in splits.values() for k,p in group.items()}
    if len(paths)!=len(by_id):raise ValueError('Record IDs must be unique across splits')
    if any(Path(k).name!=k or k in ('','.','..') for k in by_id):
        raise ValueError('Record IDs must be filenames')
    source=cfg.get('basis_training',list(splits['training']))
    if not source or any(k not in splits['training'] for k in source):raise ValueError('Basis sources must be training records')
    ranks=cfg['ranks'];channel=cfg.get('channel',0);bc=cfg['basis']
    locked=None
    if 'locked_basis' in bc:
        old=HistoryBasis.load(root/bc['locked_basis']);locked=(old.U,old.sigma)
    basis=fit_basis([(by_id[k][0],by_id[k][1]) for k in source],bc['taps'],max(ranks),args.library,
                    gap=bc.get('gap',1),scales=bc.get('scales'),locked=locked)
    prepared={}
    for k,(t,y,dy) in by_id.items():
        tm,b,db=basis.project(t,y,dy)
        first=len(t)-len(tm);truth=np.column_stack([y[first:,channel],dy[first:,channel]])
        lo,hi=cfg.get('windows',{}).get(k,[tm[0],tm[-1]])
        sl=(tm>=lo)&(tm<=hi)
        if np.count_nonzero(sl)<2:raise ValueError(f'{k} lacks a usable forecast interval')
        prepared[k]=dict(time=tm[sl],b=b[sl],db=db[sl],truth=truth[sl])
    result_dir=args.output.resolve();result_dir.mkdir(parents=True,exist_ok=False)
    basis.save(result_dir/'basis.npz')
    candidates=[dict(c) for c in cfg['candidates']]
    for c in candidates:
        if c.get('family')=='local':c['library']=str(args.library.resolve())
    study=rank_sweep(basis,ranks,*[{k:prepared[k] for k in splits[r]} for r in splits],
                     candidates,channel=channel,solver=cfg.get('solver'),validation_weights=cfg.get('validation_weights'))
    for result,field,curves in study:
        dest=result_dir/f"rank-{result['rank']}";dest.mkdir()
        if field is not None:
            if 'controls' in cfg:
                result['controls']=run_controls(cfg['controls'],result['rank'],splits,prepared,field,curves,
                    basis.endpoint(result['rank'],channel),cfg.get('solver',{}))
            field.save(dest/('model.json' if hasattr(field,'spec') else 'model.npz'))
            for label,curve in curves.items():
                np.savez(dest/(label+'.npz'),**curve)
        (dest/'scores.json').write_text(json.dumps(result,indent=2,allow_nan=False),encoding='utf-8')
    (result_dir/'study.json').write_text(json.dumps(dict(manifest=cfg,source=str(args.manifest.resolve()),
        schema='dynsyslab.miscellaneous-computes/result/v1'),indent=2),encoding='utf-8')
if __name__=='__main__':main()
