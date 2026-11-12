"""Coordinate chart, derivative pushforward and portable modal ODE with explicit units and field space."""
from pathlib import Path
import json
import numpy as np
from . import field_codec
from .tapped_history import full_tap_svd
from .models import SparseField
from .local_field import LocalField
from .compact_field import CompactField

class ModalChart:
    """Encode observation windows in signed modes and decode them with the saved lag and scaling convention."""
    def __init__(self,encoder,decoder,lags,channel_scales,working_gain=1.):
        self.encoder=np.asarray(encoder,float);self.decoder=np.asarray(decoder,float)
        self.lags=np.asarray(lags,float);self.scales=np.asarray(channel_scales,float)
        self.gain=float(working_gain)
        if (self.encoder.ndim!=2 or not self.encoder.shape[1] or self.decoder.shape!=self.encoder.T.shape
            or self.lags.ndim!=1 or self.scales.ndim!=1 or len(self.lags)*len(self.scales)!=len(self.encoder)
            or np.any(self.lags<0) or len(np.unique(self.lags))!=len(self.lags) or np.any(self.scales<=0)
            or self.gain<=0 or not all(np.isfinite(v).all() for v in (self.encoder,self.decoder,self.lags,self.scales,self.gain))
            or not np.allclose(self.decoder@self.encoder,np.eye(self.encoder.shape[1]),rtol=1e-7,atol=1e-7)):
            raise ValueError('Invalid signed-coordinate chart, lag order or unit conversion')
    @classmethod
    def fit(cls,blocks,lags,channel_rms,rank=None):
        b=full_tap_svd(blocks,channel_rms,rank=rank);max_rank=len(b['sigma'])
        rank=max_rank if rank is None else rank
        if not 1<=rank<=max_rank:raise ValueError('Invalid retained rank')
        return cls(b['encoder'][:,:rank],b['decoder'][:rank],lags,np.ones(len(channel_rms)),np.sqrt(b['windows']))
    def project(self,windows):
        w=np.asarray(windows,float)
        if w.shape[-2:]!=(len(self.lags),len(self.scales)):raise ValueError('Window dimensions disagree with chart')
        return (w*self.scales).reshape(*w.shape[:-2],-1)@self.encoder
    def decode(self,signed):
        x=np.asarray(signed)@self.decoder
        return x.reshape(*x.shape[:-1],len(self.lags),len(self.scales))/self.scales
    def record(self,time,windows,window_rates):
        from .sampling import clock
        t=clock(time);b=self.project(windows);db=self.project(window_rates)
        if b.ndim!=2 or b.shape!=db.shape or len(b)!=len(t):raise ValueError('Unaligned projected derivative labels')
        return dict(time=t,b=self.gain*b,db=self.gain*db)

class ModalModel:
    """Autonomous field and frozen chart with signed-state evaluation and pickle-free persistence."""
    schema='dynsyslab.miscellaneous-computes/modal-model/v1'
    def __init__(self,chart,field,field_space='direct',metadata=None,feature_asinh=0.):
        if field_space not in ('direct','decoded'):raise ValueError('Unknown field space')
        if field_space=='decoded' and chart.encoder.shape[0]!=chart.encoder.shape[1]:
            raise ValueError('Decoded dynamics require an invertible chart')
        if not np.isfinite(feature_asinh) or feature_asinh<0:raise ValueError('Invalid feature input scale')
        self.chart,self.field,self.space=chart,field,field_space
        self.metadata={} if metadata is None else dict(metadata)
        self.feature_asinh=float(feature_asinh)
    def rhs(self,signed):
        v=np.asarray(signed,float)
        if self.space=='decoded':return self.field(v@self.chart.decoder)@self.chart.encoder
        b=v*self.chart.gain
        inputs=np.arcsinh(b/self.feature_asinh) if self.feature_asinh else b
        return self.field(inputs)/self.chart.gain
    def save(self,path):
        a={};spec=field_codec.encode(self.field,a);c=self.chart
        with Path(path).open('xb') as stream:
            np.savez_compressed(stream,schema=self.schema,encoder=c.encoder,decoder=c.decoder,lags=c.lags,
                scales=c.scales,gain=c.gain,field_space=self.space,feature_asinh=self.feature_asinh,
                metadata=json.dumps(self.metadata,allow_nan=False),field=json.dumps(spec,allow_nan=False),**a)
    @classmethod
    def load(cls,path):
        with np.load(path,allow_pickle=False) as z:
            schema=str(z['schema'])
            if schema==cls.schema:
                c=ModalChart(z['encoder'],z['decoder'],z['lags'],z['scales'],float(z['gain']))
                return cls(c,field_codec.decode(json.loads(str(z['field'])),z),str(z['field_space']),
                           json.loads(str(z['metadata'])),float(z['feature_asinh']))
            if schema=='dynsyslab.rpe-full-delay-sindy/v1':
                from .compat_rpe import HistoricalPressureField
                scale=float(z['source_scale_per_kpa'])
                c=ModalChart(z['encoder'],z['decoder'],z['lags_tau'],[scale,scale],np.sqrt(float(z['windows'])))
                return cls(c,HistoricalPressureField(json.loads(str(z['law'])),c.lags),'decoded',
                           dict(source_schema=schema))
            if schema!='dynsyslab.direct-modal-sindy/v1':raise ValueError(f'Unsupported modal model schema {schema}')
            # This historical RPE format stored a 5 mm intrinsic pressure factor implicitly.
            c=ModalChart(z['encoder'],z['decoder'],z['lags_tau'],[5.,5.],float(z['rootN']))
            family=str(z['family']);extra=0.
            if family=='global':field=SparseField(json.loads(str(z['spec'])),z['coef'])
            elif family=='local':
                field=LocalField({k:z[k] for k in ('center','mean','P','Q','coef','scale')},float(z['decay']),int(z['neighbors']))
                extra=float(z['asinh_input'])
            elif family=='compact-linear':
                field=CompactField(z['centers'],z['coef'],z['linear'],int(z['neighbors']),float(z['asinh_input']))
            else:raise ValueError(f'Unsupported historical family {family}')
            return cls(c,field,metadata=dict(source_schema=schema,selection=json.loads(str(z['selection']))),feature_asinh=extra)
