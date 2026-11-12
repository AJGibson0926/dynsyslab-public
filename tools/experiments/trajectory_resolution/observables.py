"""Prepare independently scaled display panels for arbitrary multichannel observables."""
import numpy as np


def observable_panels(time, arrays, terms, channel_labels):
    """Yield (title, unit, curves) for (array key, title, SI unit) specifications."""
    time=np.asarray(time)
    if time.ndim!=1 or not len(time) or not np.isfinite(time).all() or np.any(np.diff(time)<=0):
        raise ValueError('Observables require a finite increasing clock')
    for key,title,unit in terms:
        values=np.asarray(arrays[key])
        if values.ndim==1:
            values=values[:,None]
        if values.shape!=(len(time),len(channel_labels)) or not np.isfinite(values).all():
            raise ValueError('Invalid observable shape or values for '+key)
        peak=float(np.max(abs(values)))
        exponent=int(np.clip(3*np.floor(np.log10(peak)/3),-12,9)) if peak and unit else 0
        prefix={-12:'p',-9:'n',-6:'µ',-3:'m',0:'',3:'k',6:'M',9:'G'}[exponent]
        label=title if peak else title+' (identically zero on this record)'
        yield label,prefix+unit,[(name,time,values[:,i]/10.**exponent) for i,name in enumerate(channel_labels)]
