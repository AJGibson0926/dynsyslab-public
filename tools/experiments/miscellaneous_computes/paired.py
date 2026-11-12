"""Phase- and amplitude-dependent sparse models for paired real modal coordinates."""
import numpy as np
from .models import records_checked
from .features import thresholded_least_squares

def features(b,j,h,phase):
    """Construct the amplitude and relative-phase features for one complex modal pair."""
    u=b[:,::2]+1j*b[:,1::2];p=abs(u)**2
    E=p.mean(axis=1)/2;u=u/np.sqrt(p+.008)
    q=p[:,j];p0=p[:,0]
    f=[np.ones_like(q),E/(E+.15),E/(E+1.5),E/(E+8),p0/(1+p0),q/(1+q)]
    if phase:
        a=u[:,0];base=a**abs(h)
        offset=u[:,j]*(base.conj() if h>=0 else base)
        f.extend([a.real,a.imag,(a*a).real,(a*a).imag,offset.real,offset.imag])
    return np.column_stack(f)

class PairedField:
    """Autonomous vector field that learns radial and angular rates for paired coordinates."""
    def __init__(self,harmonics,phase,coef):
        self.harmonics=np.asarray(harmonics,dtype=int);self.phase=bool(phase);self.coef=np.asarray(coef)
        if self.coef.shape!=(len(self.harmonics),12 if phase else 6,2):
            raise ValueError('Expected two coefficient columns per modal pair')
    def __call__(self,b):
        x=np.atleast_2d(b)
        if x.shape[1]!=2*len(self.harmonics):raise ValueError('Paired state size mismatch')
        result=np.empty_like(x,dtype=float)
        for j,h in enumerate(self.harmonics):
            a,w=(features(x,j,h,self.phase)@self.coef[j]).T
            u,v=x[:,2*j:2*j+2].T
            result[:,2*j]=a*u-w*v;result[:,2*j+1]=w*u+a*v
        return result[0] if np.ndim(b)==1 else result
    def save(self,path):
        with open(path,'xb') as stream:
            np.savez(stream,schema='dynsyslab.miscellaneous-computes/paired-field/v1',harmonics=self.harmonics,phase=self.phase,coef=self.coef)
    @classmethod
    def load(cls,path):
        with np.load(path) as z:
            if str(z['schema'])!='dynsyslab.miscellaneous-computes/paired-field/v1':raise ValueError('Unsupported paired schema')
            return cls(z['harmonics'],bool(z['phase']),z['coef'])
    @classmethod
    def fit(cls,records,*,harmonics=None,phase=True,weak_windows=(),instant_weight=1.,ridge=3e-5,threshold=.003):
        from scipy.integrate import cumulative_trapezoid
        records=records_checked(records);n=records[0]['b'].shape[1]
        if n%2 or not 0<=instant_weight<=1 or (instant_weight<1 and not weak_windows):
            raise ValueError('Supply modal pairs and valid derivative/integral weights')
        if harmonics is None:
            x=records[0]['b'];angle=np.unwrap(np.angle(x[:,::2]+1j*x[:,1::2]),axis=0)
            turns=angle[-1]-angle[0]
            if abs(turns[0])<1e-8:raise ValueError('Specify harmonics for a nonwinding reference pair')
            harmonics=np.rint(turns/turns[0]).astype(int).clip(-24,24)
        if len(harmonics)!=n//2:raise ValueError('Specify one winding number per pair')
        coeff=[]
        for j,h in enumerate(harmonics):
            gram=0.
            for r in records:
                t,b,db=r['time'],r['b'],r['db']
                f=features(b,j,h,phase);x=b[:,2*j:2*j+2];v=db[:,2*j:2*j+2]
                power=np.sum(x*x,axis=1);radius=np.sqrt(np.maximum(power,1e-5))
                angular=x[:,0]*v[:,1]-x[:,1]*v[:,0]
                data=np.column_stack([f*radius[:,None],np.sum(x*v,axis=1)/radius,angular/radius])
                gram+=instant_weight*(data.T@data)
                if instant_weight<1:
                    S=cumulative_trapezoid(f,t,axis=0,initial=0)
                    polar=np.column_stack([.5*np.log(power+1e-10),np.unwrap(np.angle(x[:,0]+1j*x[:,1]))])
                    for k in weak_windows:
                        if not isinstance(k,int) or not 1<=k<len(b):raise ValueError('Invalid weak window length')
                        X=np.column_stack([S[k:]-S[:-k],polar[k:]-polar[:-k]])/(t[k:]-t[:-k])[:,None]
                        X*=(power[k:]*power[:-k])[:,None]**.25
                        gram+=(1-instant_weight)/len(weak_windows)*(X.T@X)
            coeff.append(thresholded_least_squares(gram,2,ridge,threshold).T)
        return cls(harmonics,phase,coeff)
