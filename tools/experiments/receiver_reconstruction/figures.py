"""Render multi-ring estimates, per-ring errors and a traced rebound event."""
import argparse
import csv
from functools import partial
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator
from diagnose import read_run
from study import load, measurements, TERMS
from phase_data import DisplayVertices, interval, phase_windows, read_trace
from observables import observable_panels

BLUE,ORANGE,PURPLE,GREEN='#245f7e','#bd5626','#7955a0','#368474'


def add_chart(pages,page,axis,title,xlabel,ylabel,series,kind='time',*,select):
    traces=[]
    for name,x,y,color,style in series:
        keep=np.arange(len(x)) if len(x)<=4096 else select(x,y)
        if axis is not None:
            axis.plot(x[keep],y[keep],label=name,c=color,lw=1.65 if style=='-' else 1.45,ls=style)
        traces.append(dict(label=name,x=x[keep].tolist(),y=y[keep].tolist(),color=color,
                           source_samples=len(x),displayed_vertices=len(keep),
                           line_style='solid' if style=='-' else 'dash'))
    if axis is not None:
        axis.set(title='',xlabel=xlabel,ylabel=ylabel)
        axis.set_title(title,loc='left',fontsize=11.4,weight='bold')
        axis.grid(color='#dce5e8',lw=.65);axis.set_axisbelow(True);axis.xaxis.set_major_locator(MaxNLocator(6))
        axis.margins(x=.01)
        axis.legend(frameon=False,fontsize=8,loc='best')
    pages.append(dict(kind=kind,title=title,x_label=xlabel,y_label=ylabel,value_scale='linear',
                      dashboard_page=page,traces=traces))


def sheet(title,note):
    fig,axes=plt.subplots(2,2,figsize=(11,6.6))
    fig.subplots_adjust(left=.095,right=.979,top=.80,bottom=.105,hspace=.67,wspace=.29)
    fig.text(.068,.93,title,fontsize=19,weight='bold',color='#243b49')
    fig.text(.068,.871,note,fontsize=10.8,color='#425866')
    return fig,axes


def render(folder,diagnosis,dest,library,replays=None,detail_window=None,compare_vertices=None,term_signals=None):
    chart=partial(add_chart,select=DisplayVertices(library))
    stats,d=read_run(folder); info=json.loads(diagnosis.read_text())
    held=load(Path(stats['held_out']));train=load(Path(stats['training']))
    t=d['time'];b=d['bounds'];X=d['truth']; pages=[]
    near=stats['cases']['all']['estimates']
    N=d[near['two']['computed_as']+'_two']
    C=d[near['aligned']['computed_as']+'_aligned'] if 'scores' in near['aligned'] else None
    nshow=min(8,len(b)-1)
    if nshow<1:raise ValueError('A ring interval is required to plot this study')
    early,last,full=phase_windows(folder,held,N,b,nshow,replays)
    nrings=info['ring_intervals']; labels=['Radius r','Radial velocity dr/dτ']
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'path.simplify':False,'axes.spines.top':False,
                         'axes.spines.right':False,'axes.edgecolor':'#8d9da6'})
    fig,ax=sheet('State estimates over successive rings',
                 f'Held-out ringdown   ·   Upper panels show the first {nshow} ring intervals. Lower panels cover the full record.')
    with (folder/'by-ring.csv').open() as f:rows=list(csv.DictReader(f))
    first=(t>=b[0])&(t<b[nshow])
    for j,(name,key) in enumerate(zip(['Radius','Radial velocity'],['r_minus_1_error_pct','velocity_error_pct'])):
        chart(pages,0,ax[0,j],f'{name} over {nshow} rings','Near-receiver emission time τ',labels[j],[
            ('RPE',early.time,early.state[:,j],BLUE,'-'),
            ('Two pressures, all terms',early.query_time,early.prediction[:,j],ORANGE,'--')])
        lines=[]
        for mode,legend,col,sty in [('two','Two pressures',ORANGE,'-'),('aligned','With aligned history',PURPLE,'--')]:
            selected=[r for r in rows if r['estimator']==mode and r['case']=='all']
            if not selected:continue
            lines.append((legend,np.array([int(r['ring']) for r in selected]),
                          np.array([float(r[key]) for r in selected]),col,sty))
        chart(pages,1,ax[1,j],name+' error by ring','Ring number','Relative RMS error (%)',lines)
        for page,mask,title in [(2,t>=t[0],'Full record'),(3,(t>=b[-nshow-1])&(t<b[-1]),f'Last {nshow} rings')]:
            w=full if page==2 else last
            curves=[('RPE',w.time,w.state[:,j],BLUE,'-'),('Two pressures',w.query_time,w.prediction[:,j],ORANGE,'--')]
            if C is not None:curves.append(('With aligned history',t[mask],C[mask,j],PURPLE,'--'))
            chart(pages,page,None,title+' | '+name,'Emission time τ',labels[j],curves)
    fig.savefig(dest/'multi-ring.png',dpi=155);plt.close(fig)
    fig,ax=sheet('Why the rebound estimate jumps',
                 f'Noiseless pressure signals   ·   Each state estimate uses the {stats["k"]} nearest training samples.')
    i=info['example_index'];at=t[i]; span=b[1]-b[0]
    local=(t>=max(t[0],info['collapse_tau']+.015*span))&(t<min(b[1],at+.12*span))
    zoom=abs(t-(t[i]+t[i+1])/2)<12*(t[i+1]-t[i])
    local_ref=interval(early.time,t[local][0],t[local][-1])
    zoom_ref=interval(early.time,t[zoom][0],t[zoom][-1])
    local_query=interval(early.query_time,t[local][0],t[local][-1])
    zoom_query=interval(early.query_time,t[zoom][0],t[zoom][-1])
    raw=d['monopole_observations']/held['ranges_m']
    chart(pages,4,ax[0,0],'Receiver pressures after the first collapse','Emission time τ','Pressure (kPa)',[
        (f'{D*1000:g} mm',t[local],raw[local,j]/1000,col,'-') for j,(D,col) in enumerate(zip(held['ranges_m'],[BLUE,GREEN]))])
    chart(pages,4,ax[0,1],'Velocity during the rebound','Emission time τ','Radial velocity dr/dτ',[
        ('RPE',early.time[local_ref],early.state[local_ref,1],BLUE,'-'),
        ('Two pressures, all terms',early.query_time[local_query],early.prediction[local_query,1],ORANGE,'--')])
    chart(pages,5,ax[1,0],'Two successive readings at the jump','Emission time τ','Radial velocity dr/dτ',[
        ('RPE',early.time[zoom_ref],early.state[zoom_ref,1],BLUE,'-'),
        ('Two pressures, all terms',early.query_time[zoom_query],early.prediction[zoom_query,1],ORANGE,'--')])
    for a in [ax[0,0],ax[0,1],ax[1,0]]:a.axvspan(t[i],t[i+1],color='#577889',alpha=.15)
    for j in (i,i+1):ax[1,0].scatter(t[j],d['monopole_two'][j,1],s=26,c=ORANGE,zorder=5)
    ids=d['neighbor_rows'][zoom];weights=d['weights'][zoom]
    rings=np.searchsorted(train['ring_bounds'],d['train_time'][ids],side='right')
    top=sorted(np.unique(rings),key=lambda idx: -np.sum(weights[rings==idx]))[:3]
    group=[(f'Training ring {idx}',t[zoom],100*np.sum(weights*(rings==idx),axis=1),color,'-')
           for idx,color in zip(top,[BLUE,PURPLE,GREEN])]
    chart(pages,5,ax[1,1],'Which training rings supply the estimate?','Emission time τ','Weight (%)',group)
    ax[1,1].axvspan(t[i],t[i+1],color='#577889',alpha=.15)
    fig.savefig(dest/'rebound-variation.png',dpi=155);plt.close(fig)
    for data,role in [(train,'Training'),(held,'Held-out')]:
        raw,_=measurements(data,(0,1,2,3)); raw/=data['ranges_m']
        for j,D in enumerate(data['ranges_m']):
            chart(pages,6,None,f'{role} | {1000*D:g} mm','Near-receiver emission time τ','Acoustic pressure (kPa)',[
                ('Total pressure',data['time'],raw[:,j]/1000,BLUE,'-')])
        for j,D in enumerate(data['ranges_m']):
            chart(pages,7,None,f'{role} | {D*1000:g} mm','Near-receiver emission time τ','Correction (µPa)',[
                (label,data['time'],data[term][:,j]*1.e6,color,style)
                for term,label,color,style in zip(TERMS[1:],['Spherical near field','Dipole induction','Dipole radiation'],
                                                  [ORANGE,BLUE,PURPLE],['-','-','--'])])
    note=('RPE, pressures and estimates use every solver output in these windows.' if replays else
          'RPE uses every solver output. Pressures and estimates use the receiver clock.')
    fig,ax=sheet('State and pressure phase portraits','Held-out ringdown   ·   '+note)
    late=(t>=b[-nshow-1])&(t<b[-1])
    pressure_labels=[f'{D*1000:g} mm pressure (kPa)' for D in held['ranges_m']]
    for column,(mask,label) in enumerate([(first,f'First {nshow} ring intervals'),(late,f'Last {nshow} rings'),
                                          (np.ones(len(t),bool),'Full ringdown')]):
        state_axis=ax[0,column] if column<2 else None
        pressure_axis=ax[1,column] if column<2 else None
        w=(early,last,full)[column]
        chart(pages,8 if column<2 else 10,state_axis,label+' | bubble state','Radius r','Radial velocity dr/dτ',[
            ('True RPE state',w.state[:,0],w.state[:,1],BLUE,'-'),
            ('Two-pressure estimate',w.prediction[:,0],w.prediction[:,1],ORANGE,'--')],kind='phase')
        chart(pages,9 if column<2 else 10,pressure_axis,label+' | receiver pressures',*pressure_labels,[
            ('Monopole',w.leading[:,0],w.leading[:,1],BLUE,'-'),
            ('All terms',w.total[:,0],w.total[:,1],GREEN,'--')],kind='phase')
        if C is not None:
            chart(pages,11,None,label+' | aligned pressure history','Radius r','Radial velocity dr/dτ',[
                ('True RPE state',w.state[:,0],w.state[:,1],BLUE,'-'),
                ('With aligned history',C[mask,0],C[mask,1],PURPLE,'--')],kind='phase')
    fig.savefig(dest/'phase-space.png',dpi=155);plt.close(fig)
    names=[f'First {nshow} rings',f'Errors by ring · {nrings} intervals','Full state record',f'Last {nshow} rings',
           'Rebound | input and output','Rebound | neighbor weights','Both pressure records','Each correction',
           'State phase | early and late','Pressure plane | early and late','Full phase portraits']
    if C is not None:names.append('State phase | aligned history')
    ablation=len(names);names.append('State phase | individual corrections')
    for mask,section,w in [(first,'Early rings',early),(late,'Late rings',last)]:
        curves=[('RPE',w.state[:,0],w.state[:,1],BLUE,'-')]
        for case,legend,col,sty in [('monopole','Monopole','#546671','-'),('near_field','Add spherical',ORANGE,'--'),
                                    ('induction','Add dipole induction',PURPLE,'--'),('radiation','Add dipole radiation',GREEN,'--')]:
            key=stats['cases'][case]['estimates']['two']['computed_as']+'_two'
            curves.append((legend,d[key][mask,0],d[key][mask,1],col,sty))
        chart(pages,ablation,None,section+' | corrections','Radius r','Radial velocity dr/dτ',curves,kind='phase')
    reference_page=len(names);names.append('Both ringdowns | physical and receiver phase')
    for source,role in [(train,'Training'),(held,'Held-out')]:
        if role=='Held-out':
            state=full.state
        else:
            native_t,state=read_trace(source['meta']['source'])
            state=state[interval(native_t,source['time'][0],source['time'][-1])]
        chart(pages,reference_page,None,role+' | full physical state','Radius r','Radial velocity dr/dτ',[
            ('RPE',state[:,0],state[:,1],BLUE,'-')],kind='phase')
        p,_=measurements(source,range(len(TERMS)));p/=source['ranges_m']*1000.
        chart(pages,reference_page,None,role+' | full receiver pair',*pressure_labels,[
            ('All terms',p[:,0],p[:,1],GREEN,'-')],kind='phase')
    priority=[8,9,10]
    if term_signals:
        titles=['Leading monopole','Spherical near-field correction','Dipole induction','Dipole radiation']
        terms=list(zip(TERMS,titles,['Pa']*len(TERMS)))
        for path,base,role in zip(term_signals,(held,train),('Held-out','Training')):
            values=load(path)
            source=values['meta']
            if source['source']!=base['meta']['source'] or source['receiver_config']['stride']!=1:
                raise ValueError('Term figures require full-resolution pressure records from this study')
            for key in ('positions_m','sound_speed_m_s'):
                if source['receiver_config'][key]!=base['meta']['receiver_config'][key]:
                    raise ValueError('Term figure receiver configuration differs from study')
            page=len(names);names.append(role+' | four acoustic observables')
            priority.append(page)
            fig,axes=sheet(role+' ringdown | acoustic pressure terms',
                'First collapse and rebound   ·   Each receiver uses its retarded emission. Each panel has its own pressure scale.')
            channels=[f'{D*1000:g} mm' for D in values['ranges_m']]
            for ax,(label,unit,curves) in zip(axes.flat,observable_panels(values['time'],values,terms,channels)):
                series=[(*curve,color,sty) for curve,color,sty in zip(curves,(BLUE,ORANGE),('-','--'))]
                chart(pages,page,ax,label,'Near-receiver emission time τ','Pressure ('+unit+')',series)
                if not any(np.any(curve[2]) for curve in curves):ax.set_yticks([0])
            fig.savefig(dest/(role.lower()+'-observables.png'),dpi=155);plt.close(fig)
    if detail_window:
        if detail_window[1]<=detail_window[0] or detail_window[0]<full.time[0] or detail_window[1]>full.time[-1]:
            raise ValueError('Detail window must lie inside the observed record')
        idx=interval(full.time,*detail_window)
        if len(full.time[idx])<2:raise ValueError('Detail window requires at least two solver outputs')
        detail=len(names);names.append('Resolution check | phase and time')
        priority.insert(1,detail)
        nt,xx=full.time[idx],full.state[idx]
        phase=[('Native RPE',xx[:,0],xx[:,1],BLUE,'-')]
        velocity=[('Native RPE',nt,xx[:,1],BLUE,'-')]
        if compare_vertices:
            st,ss=read_trace(compare_vertices,state_key='state')
            if st[0]>nt[0] or st[-1]<nt[-1]:raise ValueError('Comparison vertices must span the detail window')
            st1=np.r_[nt[0],st[(st>nt[0])&(st<nt[-1])],nt[-1]]
            r,v=(np.interp(st1,st,component) for component in ss.T)
            phase.append(('Previous plot',r,v,ORANGE,'--'))
            velocity.append(('Previous plot',st1,v,ORANGE,'--'))
        chart(pages,detail,None,'Sampling through the rebound | phase','Radius r','Radial velocity dr/dτ',phase,kind='phase')
        chart(pages,detail,None,'Sampling through the rebound | velocity','Emission time τ','Radial velocity dr/dτ',velocity)
    order=priority+[i for i in range(len(names)) if i not in priority]
    for f in pages:f['dashboard_page']=order.index(f['dashboard_page'])
    names=[names[i] for i in order]
    payload=dict(schema='dynsyslab.viewer.payload/v1',source=str((folder/'results.json').resolve()),view='phase',
                 tile=dict(name='Two receivers | correction study',columns=2),
                 page_labels=names,
                 figures=sorted(pages,key=lambda f:f['dashboard_page']))
    (dest/'payload.json').write_text(json.dumps(payload,separators=(',',':')),encoding='utf-8')
    (dest/'catalog.json').write_text(json.dumps(dict(schema='dynsyslab.viewer.catalog/v1',items=[
        dict(id='receiver-corrections',label='Receiver corrections and rebound diagnostics',kind='trajectory',source='payload.json',views=['phase'])])))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('study','diagnosis','out','library'):parser.add_argument('--'+key,type=Path,required=True)
    parser.add_argument('--replays',nargs=2,type=Path,metavar=('EARLY','LATE'))
    parser.add_argument('--detail-window',nargs=2,type=float,metavar=('START','STOP'))
    parser.add_argument('--compare-vertices',type=Path)
    parser.add_argument('--term-signals',nargs=2,type=Path,metavar=('HELD','TRAINING'))
    args=parser.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    render(args.study,args.diagnosis,args.out,args.library,args.replays,args.detail_window,args.compare_vertices,args.term_signals)
