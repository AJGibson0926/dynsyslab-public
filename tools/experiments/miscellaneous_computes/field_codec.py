"""Versioned, pickle-free numeric fields for deployable side-car experiments."""
import numpy as np
from .models import SparseField
from .paired import PairedField
from .compact_field import CompactField
from .local_field import LocalField
from .structured import SumMap,SecondOrderField,PairNetwork,InteriorField,EndpointAcceleration

def encode(field,arrays):
    """Record field structure in a JSON tree and append its numeric parameters to the supplied array mapping."""
    from .compat_rpe import HistoricalPressureField
    def store(value):
        name=f'array_{len(arrays)}';arrays[name]=np.asarray(value);return name
    child=lambda f:encode(f,arrays)
    if isinstance(field,EndpointAcceleration):
        m=field.model;c=m.chart
        return dict(type='endpoint-model',encoder=store(c.encoder),decoder=store(c.decoder),lags=store(c.lags),
                    scales=store(c.scales),gain=c.gain,space=m.space,asinh_scale=m.feature_asinh,
                    metadata=m.metadata,field=child(m.field))
    if isinstance(field,HistoricalPressureField):return dict(type='rpe-legacy',law=field.law,lags=store(field.lags))
    if isinstance(field,SparseField):return dict(type='sparse',spec=field.spec,C=store(field.coef))
    if isinstance(field,CompactField):
        return dict(type='compact',centers=store(field.features.centers),C=store(field.coefficients),
                    linear=store(field.linear),k=field.features.k,scale=field.features.scale)
    if isinstance(field,LocalField):
        return dict(type='local',arrays={k:store(v) for k,v in field.arrays.items()},k=field.neighbors,decay=field.decay)
    if isinstance(field,PairedField):
        return dict(type='oscillator',harmonics=store(field.harmonics),phase=field.phase,C=store(field.coef))
    if isinstance(field,SumMap):return dict(type='sum',fields=[child(f) for f in field.fields])
    if isinstance(field,SecondOrderField):return dict(type='second-order',acceleration=child(field.acceleration))
    if isinstance(field,PairNetwork):
        return dict(type='network',lags=store(field.lags),weights=store(field.weights),pairs=field.pairs,fields=[child(f) for f in field.fields])
    if isinstance(field,InteriorField):
        return dict(type='interior',outer=child(field.outer),drive=child(field.drive),gamma=store(field.gamma),
                    feedback=None if field.feedback is None else child(field.feedback))
    raise ValueError(f'Unsupported portable field {type(field).__name__}')

def decode(spec,arrays):
    """Rebuild the frozen evaluator, including embedded endpoint models, from its operator tree and numeric arrays."""
    typ=spec['type'];a=lambda key:arrays[spec[key]];child=lambda s:decode(s,arrays)
    if typ=='endpoint-model':
        from .modal import ModalChart,ModalModel
        c=ModalChart(a('encoder'),a('decoder'),a('lags'),a('scales'),spec['gain'])
        return EndpointAcceleration(ModalModel(c,child(spec['field']),spec['space'],spec['metadata'],spec['asinh_scale']))
    if typ=='rpe-legacy':
        from .compat_rpe import HistoricalPressureField
        return HistoricalPressureField(spec['law'],a('lags'))
    if typ=='sparse':return SparseField(spec['spec'],a('C'))
    if typ=='compact':return CompactField(a('centers'),a('C'),a('linear'),spec['k'],spec['scale'])
    if typ=='local':return LocalField({k:arrays[v] for k,v in spec['arrays'].items()},spec['decay'],spec['k'])
    if typ=='oscillator':return PairedField(a('harmonics'),spec['phase'],a('C'))
    if typ=='sum':return SumMap([child(f) for f in spec['fields']])
    if typ=='second-order':return SecondOrderField(child(spec['acceleration']))
    if typ=='network':return PairNetwork(a('lags'),a('weights'),spec['pairs'],[child(f) for f in spec['fields']])
    if typ=='interior':return InteriorField(child(spec['outer']),child(spec['drive']),a('gamma'),
                                            None if spec['feedback'] is None else child(spec['feedback']))
    raise ValueError(f'Unsupported portable field type {typ}')
