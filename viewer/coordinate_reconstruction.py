"""Presentation-only figures for supplied coordinate decoders and error studies."""
from __future__ import annotations


def reconstruction_figure(records, *, time_limit=None, labels=('x₁','x₂')):
    """Each row: time, reference Nx2, decoded Nx2, label; optional ode Nx2."""
    import numpy as np
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(len(records),3,figsize=(16,3.5*len(records)),squeeze=False)
    for row,(record,axs) in enumerate(zip(records,axes)):
        t=np.asarray(record['time']);ref=np.asarray(record['reference']);decoded=np.asarray(record['decoded'])
        keep=t<=time_limit if time_limit is not None else np.ones(len(t),dtype=bool)
        for c in range(2):
            if 'ode' in record:axs[c].plot(t[keep],np.asarray(record['ode'])[keep,c],color='#b57b00',lw=3,label='Independent ODE')
            axs[c].plot(t[keep],ref[keep,c],color='#c52b24',lw=2,label='Closed form')
            axs[c].plot(t[keep],decoded[keep,c],color='#0072b2',lw=1.2,ls='--',label='Character reconstruction')
            axs[c].set(xlabel='Time',ylabel=labels[c])
        axs[2].plot(ref[:,0],ref[:,1],color='#c52b24',lw=2.8,label='Closed form')
        axs[2].plot(decoded[:,0],decoded[:,1],color='#0072b2',lw=1.2,ls='--',label='Character reconstruction')
        axs[2].set(xlabel=labels[0],ylabel=labels[1])
        for ax in axs:ax.grid(alpha=.18);ax.set_title(record['label'])
    handles,text=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,text,loc='upper center',ncol=3)
    fig.tight_layout(rect=(0,0,1,.955))
    return fig


def character_mass_figure(records, *, title='Character contributions across finite cycles'):
    """Plot native conductor_masses; blank cells denote absent characters.

    Denominator one (the constant character) is omitted. Entries are supplied
    vector-observable squared masses, not singular values or stability exponents.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    orders=sorted({int(g['denominator']) for r in records for g in r['conductor_masses'] if g['denominator']>1})
    values=np.zeros((len(orders),len(records)))
    for j,r in enumerate(records):
        for group in r['conductor_masses']:
            if group['denominator']>1:values[orders.index(group['denominator']),j]=group['mass']
    positive=values[values>0]
    if not len(positive):raise ValueError('No nonconstant character mass supplied')
    fig,ax=plt.subplots(figsize=(11,6))
    cm=plt.get_cmap('viridis').copy();cm.set_bad('#f0f0f0')
    im=ax.imshow(np.ma.masked_where(values<=0,values),origin='lower',aspect='auto',cmap=cm,
                 norm=LogNorm(vmin=positive.min(),vmax=positive.max()))
    for i in range(len(orders)):
        for j in range(len(records)):
            rgb=cm(im.norm(values[i,j]))[:3] if values[i,j]>0 else (1,1,1)
            linear=[c/12.92 if c<=.04045 else ((c+.055)/1.055)**2.4 for c in rgb]
            luminance=sum(c*w for c,w in zip(linear,(.2126,.7152,.0722)))
            ax.text(j,i,f'{values[i,j]:.1e}' if values[i,j]>0 else 'absent',ha='center',va='center',
                    color='white' if luminance<.179 else 'black',fontsize=9)
    ax.set(xticks=range(len(records)),xticklabels=[r['label'] for r in records],
           yticks=range(len(orders)),yticklabels=[str(o) for o in orders],
           ylabel='Character order / reduced denominator',xlabel='Physical cycle and parameter',title=title)
    fig.colorbar(im,ax=ax,label='Squared coefficient mass (constant part omitted)')
    fig.tight_layout();return fig


def error_comparison_figure(x, curves, *, x_label='Resolution', y_label='Error (stated norm)', title='Reconstruction diagnostics'):
    """Draw supplied positive errors, with explicit series names and no fitting."""
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(10,5))
    colors=['#c52b24','#0072b2','#b57b00','#00865a']
    for index,(label,values) in enumerate(curves.items()):
        color=colors[index%len(colors)]
        points=[(a,b) for a,b in zip(x,values) if b>0]
        if points:ax.loglog(*zip(*points),'-o',label=label,color=color)
    ax.set(xlabel=x_label,ylabel=y_label,title=title)
    ax.legend();ax.grid(which='both',alpha=.2);fig.tight_layout();return fig
