"""Standalone CLI for modal studies, replay, scoring, reconstruction and native Viewer output."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import numpy as np

if __package__ in ('',None):
    sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from miscellaneous_computes.manifest import publish_directory,split_paths,file_id
from miscellaneous_computes.modal import ModalModel
from miscellaneous_computes.modal_study import run_manifest
from miscellaneous_computes.replay import replay
from miscellaneous_computes.diagnostics import forecast
from miscellaneous_computes.equations import modal_equations
from miscellaneous_computes.crossings import candidates,refine
from miscellaneous_computes.reconstruction import compare_reconstructions
from miscellaneous_computes.comparison import comparison_pages
from miscellaneous_computes.viewer import dashboard,render

def read_json(path):return json.loads(Path(path).read_text(encoding='utf-8'))
def save_json(path,result):
    """Publish one UTF-8 JSON result at a new path, rejecting nonfinite values."""
    text=json.dumps(result,indent=2,allow_nan=False)
    with Path(path).open('x',encoding='utf-8') as stream:stream.write(text+'\n')

def run(args):
    """Execute a parsed standalone workflow with its declared input and output contracts."""
    if args.command=='study':
        return run_manifest(args.manifest,args.output,args.library)
    if args.command=='closure':
        query=read_json(args.query)
        if query.get('schema') not in ('dynsyslab.observable-closure-query/v1','dynsyslab.observable-closure-query/v2'):
            raise ValueError('Expected a DynSysLab observable-closure query')
        result=subprocess.run([str(args.cli.resolve()),'observable-closure',str(args.query.resolve())],
                              capture_output=True,text=True,encoding='utf-8')
        if result.returncode:
            message=result.stderr.strip() or result.stdout.strip()
            detail=message.splitlines()[0] if message else 'No diagnostic output.'
            raise RuntimeError(f'Closure query failed with {args.cli} (exit {result.returncode}). '
                               f'Use a native build supporting observable-closure. {detail}')
        save_json(args.output,json.loads(result.stdout))
        return
    if args.command=='replay':
        model=ModalModel.load(args.model)
        with np.load(args.query,allow_pickle=False) as z:
            opts=dict(method=args.method,library=args.library)
            if args.step is not None:opts['max_step']=args.step
            for key in ('initial','window'):
                if key in z:opts[key]=z[key]
            r=replay(model,z['time'],**opts)
        windows=r.pop('values')
        zero=np.flatnonzero(model.chart.lags==0)
        if len(zero)!=1:raise ValueError('Select a chart with a current measurement for CLI replay')
        with args.output.open('xb') as f:
            np.savez_compressed(f,schema='dynsyslab.miscellaneous-computes/replay/v1',model=str(args.model.resolve()),
                                windows=windows,values=windows[:,zero[0]],**r)
    elif args.command=='equations':
        save_json(args.output,modal_equations(ModalModel.load(args.model)))
    elif args.command=='score':
        with np.load(args.reference) as x,np.load(args.prediction) as y:
            a=x[args.reference_key];b=y[args.prediction_key]
            if args.reference_columns is not None:a=a[:,args.reference_columns]
            if args.prediction_columns is not None:b=b[:,args.prediction_columns]
            result=forecast(x[args.reference_time],a,y[args.prediction_time],b,
                            prominence=args.prominence,channel=args.channel,time_weighted=args.time_weighted)
        save_json(args.output,dict(schema='dynsyslab.miscellaneous-computes/forecast-score/v1',**result))
    elif args.command=='crossings':
        from scipy.interpolate import CubicSpline
        with np.load(args.source) as f:t=f['time'];x=f[args.key][:,args.axes]
        hits=candidates(t,x,args.library,min_separation=args.min_separation,limit=args.limit)
        if args.refine_with:
            with np.load(args.refine_with) as f:
                curve=CubicSpline(f['time'],f[args.key][:,args.axes],extrapolate=False)
            for h in hits:h['refinement']=refine(h,curve,velocity=curve.derivative(),tolerance=args.tolerance)
        save_json(args.output,dict(schema='dynsyslab.miscellaneous-computes/crossings/v1',source=str(args.source.resolve()),axes=args.axes,candidates=hits))
    elif args.command=='reconstruct':
        cfg=read_json(args.manifest)
        if cfg.get('schema')!='dynsyslab.miscellaneous-computes/reconstruction/v1':raise ValueError('Unsupported reconstruction study')
        root=args.manifest.resolve().parent;data={};seen=set()
        for role in ('training','test'):
            data[role]={}
            for tag,name in cfg[role].items():
                file_id(tag);p=(root/name).resolve()
                if p in seen:raise ValueError('Repeated reconstruction source')
                seen.add(p)
                with np.load(p) as z:data[role][tag]=(z['time'],z['values'],z['state'])
        result=compare_reconstructions(**data,designs=cfg['designs'],**cfg.get('options',{}))
        with publish_directory(args.output) as d:
            for design,entry in result.items():
                folder=d/file_id(design);folder.mkdir()
                for tag,case in entry['cases'].items():
                    np.savez_compressed(folder/(tag+'.npz'),**case.pop('curves'))
            save_json(d/'results.json',dict(schema='dynsyslab.miscellaneous-computes/reconstruction-results/v1',manifest=cfg,designs=result))
    elif args.command=='plot':
        cfg=read_json(args.manifest);root=args.manifest.resolve().parent
        if cfg.get('schema')!='dynsyslab.miscellaneous-computes/comparison/v1':raise ValueError('Unsupported plot manifest')
        cases={}
        for role,curves in cfg['cases'].items():
            cases[role]={}
            for name,s in curves.items():
                with np.load(root/s['file']) as z:cases[role][name]=(z[s.get('time','time')],z[s.get('values','values')])
        pages=comparison_pages(cases,cfg['models'],**cfg.get('options',{}))
        # Viewer caches can contain absolute paths; generate them at their final location.
        args.output.mkdir(parents=True,exist_ok=False)
        payload=dashboard(args.output,pages,cfg.get('title','Forecast comparison'),args.manifest.resolve(),args.library)
        for i in range(len(pages)):render(args.viewer,payload,i,args.output/f'page-{i+1}.png')

def parser():
    """Build the CLI for fitting, deployment, diagnostics, reconstruction and Viewer comparisons."""
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    for name,attrs in [('study',['manifest','output']),('replay',['model','query','output']),('equations',['model','output']),
                       ('score',['reference','prediction','output']),('crossings',['source','output']),
                       ('reconstruct',['manifest','output']),('plot',['manifest','output']),('closure',['query','output'])]:
        a=sub.add_parser(name)
        for arg in attrs:a.add_argument(arg,type=Path)
        if name in ('study','replay','crossings','plot'):a.add_argument('--library',type=Path,required=name in ('crossings','plot'))
        if name=='replay':
            a.add_argument('--method',choices=('adaptive','exact-affine','rk4'),default='adaptive');a.add_argument('--step',type=float)
        if name=='score':
            for role in ('reference','prediction'):
                a.add_argument('--'+role+'-time',default='time');a.add_argument('--'+role+'-key',default='values')
                a.add_argument('--'+role+'-columns',nargs='+',type=int)
            a.add_argument('--prominence',type=float);a.add_argument('--channel',type=int,default=0);a.add_argument('--time-weighted',action='store_true')
        if name=='crossings':
            a.add_argument('--axes',nargs=2,type=int,default=[0,1]);a.add_argument('--key',default='values')
            a.add_argument('--refine-with',type=Path);a.add_argument('--tolerance',type=float,default=1e-9)
            a.add_argument('--min-separation',type=float,default=0.);a.add_argument('--limit',type=int,default=100000)
        if name=='plot':a.add_argument('--viewer',type=Path,required=True)
        if name=='closure':a.add_argument('--cli',type=Path,required=True,help='Native DynSysLab CLI with observable-closure support')
    return p

if __name__=='__main__':run(parser().parse_args())
