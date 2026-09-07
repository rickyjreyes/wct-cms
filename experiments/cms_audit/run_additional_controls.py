"""Independent unbinned/profile, width, and signal-absorption checks."""
import json,pathlib
import numpy as np
from scipy.special import logsumexp
from cms_wct.background_families import fit_background_family
from cms_wct.signature import scan_omegas,weighted_linear_sinusoid
from models import *
from run_spectrum_attacks import write


def unbinned(mass,degree,order):
    x=np.log(mass[(mass>=2)&(mass<=120)&eligible(mass)])
    nodes,weights=np.polynomial.legendre.leggauss(order);q=[];w=[]
    for lo,hi in [(2,2.9),(3.3,3.55),(3.85,8.5),(11.5,80),(100,120)]:
        lo,hi=np.log([lo,hi]);q.extend((lo+hi)/2+(hi-lo)/2*nodes);w.extend(weights*(hi-lo)/2)
    q=np.array(q);w=np.array(w)
    design=lambda z:np.polynomial.chebyshev.chebvander(2*(z-np.log(2))/np.log(60)-1,degree)[:,1:]
    def fit(X,Q,initial=None):
        b=np.zeros(Q.shape[1]) if initial is None else initial.copy();target=X.sum(0);n=len(X)
        def objective(b):return n*logsumexp(Q@b,b=w)-target@b
        for _ in range(80):
            p=w*np.exp(Q@b-logsumexp(Q@b,b=w));mean=p@Q;C=(Q-mean).T@(p[:,None]*(Q-mean));grad=n*mean-target
            step=np.linalg.lstsq(n*C,grad,rcond=1e-11)[0];scale=1.;old=objective(b)
            while objective(b-scale*step)>old+1e-8 and scale>2**-25:scale*=.5
            b-=scale*step
            if np.max(abs(Q@(scale*step)))<1e-8:break
        return b,objective(b),np.linalg.pinv(n*C,rcond=1e-11)
    X=design(x);Q=design(q);bn,ln,_=fit(X,Q);ba,la,c=fit(np.column_stack([X,np.cos(OMEGA*x-PHASE)]),np.column_stack([Q,np.cos(OMEGA*q-PHASE)]),np.r_[bn,0.])
    return dict(degree=degree,order=order,n=len(x),amplitude=ba[-1],se=np.sqrt(c[-1,-1]),delta2logL=2*(ln-la),positive_score=max(0,2*(ln-la)) if ba[-1]>0 else 0.)


def main():
    root=pathlib.Path('audit/2026-09-06');out=root/'controls';out.mkdir(exist_ok=True);ub=[];adaptive=[];inject=[]
    for label in ('H2','G1','G2'):
        mass=np.load(root/'baseline'/(label+'_masses.npz'))['masses']
        for d in (7,12,16):
            for order in (80,160):ub.append(dict(label=label,**unbinned(mass,d,order)))
        for n in (175,350,500):
            edges=np.quantile(mass[(mass>=2)&(mass<=120)],np.linspace(0,1,n+1));m=np.sqrt(edges[1:]*edges[:-1]);y=np.histogram(mass,edges)[0];g=eligible(m)
            for d in (7,12):
                f=fit_model(f'pois_poly{d}',m,y,np.diff(edges));adaptive.append(dict(label=label,bins=n,degree=d,**locked(np.log(m[g]),((y-f.mean)/np.sqrt(f.mean))[g])))
    means=np.load(root/'attacks/background_means.npz');s=np.genfromtxt(root/'baseline/H2/spectrum.csv',delimiter=',',names=True);m0=s['mass_center'];width0=s['mass_high']-s['mass_low']
    configs=[dict(name='nominal',amplitude=1.,phase=PHASE,omega=OMEGA,scale=1.,bins=350,generator='spline_s2')]
    for key,vals in [('amplitude',[.25,.5]),('phase',[PHASE+np.pi/2,PHASE+np.pi]),('omega',[OMEGA*.8,OMEGA*1.2]),('scale',[.25,4.]),('bins',[175,700]),('generator',['pois_poly12'])]:
        for v in vals:configs.append({**configs[0],key:v,'name':f'{key}_{v}'})
    for ci,c in enumerate(configs):
        edges=np.geomspace(2,120,c['bins']+1);m=np.sqrt(edges[1:]*edges[:-1]);g=eligible(m);x=np.log(m[g]);b=np.exp(np.interp(np.log(m),np.log(m0),np.log(means['H2__'+c['generator']]/width0)))*np.diff(edges)*c['scale']
        for amp in (0,c['amplitude']):
            mu=b.copy();mu[g]+=amp*np.sqrt(b[g])*np.cos(c['omega']*x-c['phase']);rng=np.random.default_rng(20260910+ci*100+int(amp>0))
            for trial in range(100):
                y=rng.poisson(mu).astype(float)
                for degree in (7,12):
                    bg,_=fit_background_family(m,y,family='chebyshev',degree=degree,excluded_windows=MASKS);r=((y-bg)/np.sqrt(bg))[g];_,best=scan_omegas(x,r,np.linspace(3.1,80,1000));fixed=weighted_linear_sinusoid(x,r,c['omega'])
                    inject.append(dict(**c,injected=amp,trial=trial,degree=degree,recovered_omega=best.omega,recovered_amplitude=fixed.amplitude,search_score=best.delta_chi2))
        print(c['name'],flush=True)
    write(out/'unbinned.csv',ub);write(out/'adaptive_width_corrected.csv',adaptive);write(out/'recovery_trials.csv',inject)
    summary=[]
    for c in configs:
        for degree in (7,12):
            rows=[r for r in inject if r['name']==c['name'] and r['degree']==degree];null=[r for r in rows if not r['injected']];sig=[r for r in rows if r['injected']];threshold=np.quantile([r['search_score'] for r in null],.95);amp=np.array([r['recovered_amplitude'] for r in sig]);freq=np.array([r['recovered_omega'] for r in sig]);power=np.mean([r['search_score']>threshold for r in sig])
            summary.append(dict(**c,degree=degree,mean_amplitude=amp.mean(),bias=amp.mean()-c['amplitude'],sd=amp.std(ddof=1),frequency_rmse=np.sqrt(np.mean((freq-c['omega'])**2)),frequency_within_0p5=np.mean(abs(freq-c['omega'])<.5),power=power,false_negative_rate=1-power,matched_null_q95=threshold))
    write(out/'recovery_summary.csv',summary)
if __name__=='__main__':main()
