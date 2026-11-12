"""Plot a dense reference and sampled vertices in phase space and on the clock."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
from resolution import compare_clocks, read_trace


def plot(reference,sampled,window,axes,labels,time_label,destination,title,sample_label='Sampled trace'):
    rt,ref=reference
    st,coarse=sampled
    report=compare_clocks(rt,ref,st,coarse,window)
    lo,hi=report['compared_time']
    native=slice(rt.searchsorted(lo),rt.searchsorted(hi,side='right'))
    # Draw the incoming and outgoing segments, including where a window cuts through one.
    interior=(st>lo)&(st<hi)
    clock=np.r_[lo,st[interior],hi]
    coords=np.column_stack([np.interp(clock,st,c) for c in coarse.T])
    markers=(st>=lo)&(st<=hi)
    a,b=axes
    with plt.rc_context({'path.simplify':False,'font.size':11,'axes.spines.right':False,'axes.spines.top':False}):
        fig,panels=plt.subplots(1,2,figsize=(11,4.7))
        fig.subplots_adjust(left=.085,right=.98,bottom=.17,top=.79,wspace=.26)
        fig.suptitle(title,x=.055,ha='left',fontsize=17,weight='bold',color='#243b49')
        for ax,nx,ox,mx,label in [(panels[0],ref[native,a],coords[:,a],coarse[markers,a],labels[0]),
                                   (panels[1],rt[native],clock,st[markers],time_label)]:
            ax.plot(nx,ref[native,b],color='#245f7e',lw=2,label=f"Native reference ({report['reference_count']:,} samples)")
            ax.plot(ox,coords[:,b],color='#ad552a',ls='--',lw=1.4,label=f"{sample_label} ({report['sample_count']:,} vertices)")
            ax.scatter(mx,coarse[markers,b],s=18,facecolor='white',edgecolor='#ad552a',zorder=3)
            ax.set(xlabel=label,ylabel=labels[1])
            ax.grid(alpha=.23)
            ax.xaxis.set_major_locator(MaxNLocator(5))
        panels[0].legend(frameon=False,fontsize=9)
        panels[1].set_xlim(lo,hi)
        fig.savefig(destination,dpi=160)
        plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('reference','sample','out'):
        p.add_argument('--'+key,type=Path,required=True)
    for key in ('reference','sample'):
        p.add_argument('--'+key+'-time', default='time')
        p.add_argument('--'+key+'-values', default='values')
    p.add_argument('--window',nargs=2,type=float,required=True)
    p.add_argument('--axes',nargs=2,type=int,default=(0,1))
    p.add_argument('--labels',nargs=2,default=('Coordinate 1','Coordinate 2'))
    p.add_argument('--time-label',default='Time')
    p.add_argument('--title',default='Sampling resolution in phase space and time')
    p.add_argument('--sample-label',default='Sampled trace')
    a=p.parse_args()
    plot(read_trace(a.reference,a.reference_time,a.reference_values),
         read_trace(a.sample,a.sample_time,a.sample_values),a.window,a.axes,a.labels,a.time_label,a.out,a.title,a.sample_label)
