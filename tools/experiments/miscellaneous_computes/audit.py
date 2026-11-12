"""Tolerance, startup-history and branch-separation controls for a proposed autonomous law."""
import numpy as np
from .study import evaluate,relative_rms

def replay_controls(field,record,rank,endpoint,*,restart_time,compare_from=None,rtol=1e-8,atol=1e-10,**solver):
    """Compare a baseline rollout to a tighter solve and a later single-state restart on common times."""
    first,original=evaluate(field,record,rank,endpoint,rtol=rtol,atol=atol,**solver)
    tight,refined=evaluate(field,record,rank,endpoint,rtol=rtol/10,atol=atol/10,**solver)
    t=record['time'];i=np.searchsorted(t,restart_time)
    if i>=len(t)-1:raise ValueError('Restart falls outside the usable trajectory')
    later={k:np.asarray(record[k])[i:] for k in ('time','b','db','truth')}
    rest,second=evaluate(field,later,rank,endpoint,rtol=rtol,atol=atol,**solver)
    report=dict(original=first,tightened=tight,restarted=rest,restart_at=float(t[i]))
    if first['complete'] and tight['complete']:
        diff=refined['predicted']-original['predicted']
        report['tolerance_change_in_reference_rms']=relative_rms(record['truth'],record['truth']+diff)
    if first['complete'] and rest['complete']:
        j=np.searchsorted(t,max(t[i],t[i] if compare_from is None else compare_from))
        if j==len(t):raise ValueError('Comparison interval is empty')
        report['common_start']=float(t[j])
        report['original_on_common']=relative_rms(record['truth'][j:],original['predicted'][j:])
        report['restart_on_common']=relative_rms(record['truth'][j:],second['predicted'][j-i:])
    return report
