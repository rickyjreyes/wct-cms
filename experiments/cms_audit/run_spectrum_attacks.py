"""Deterministic retrospective spectrum attacks; every requested grid row retained."""
import argparse,csv,json,pathlib,dataclasses
import numpy as np
from scipy.stats import chi2
from cms_wct.background_cv import frozen_background_candidates, BackgroundCandidate
from cms_wct.background_kill import PairKillConfig,evaluate_candidate,pair_cv_scores
from cms_wct.signature import scan_omegas,weighted_linear_sinusoid
from cms_wct.background import make_histogram
from models import *


def write(path, rows):
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with open(path,'w') as f:
        w=csv.DictWriter(f,keys);w.writeheader();w.writerows(rows)


def main():
    p=argparse.ArgumentParser();p.add_argument('--baseline',default='audit/2026-09-06/baseline');p.add_argument('--out',default='audit/2026-09-06/attacks');a=p.parse_args()
    out=pathlib.Path(a.out);out.mkdir(parents=True,exist_ok=True);base=pathlib.Path(a.baseline)
    specs={l:np.genfromtxt(base/l/'spectrum.csv',delimiter=',',names=True) for l in ['H1','H2','G1','G2']}
    m=specs['H2']['mass_center'];width=specs['H2']['mass_high']-specs['H2']['mass_low'];good=eligible(m);x=np.log(m)
    config=PairKillConfig(OMEGA,PHASE,MASKS)
    historical=frozen_background_candidates()+[BackgroundCandidate(f'cheb_d{d}','chebyshev',degree=d) for d in (2,3,14,16,20)]+[BackgroundCandidate(f'spline_s{s}','spline',spline_smoothing_factor=s) for s in (.25,1.5,3.,4.)]
    names=candidates();rows=[];profiles=[];means={}
    for label,s in specs.items():
        y=s['count']; means[label]={}
        for c in historical:
            try:
                r,mu=evaluate_candidate(m,y,c,config);rr=(y-mu)/np.sqrt(mu);means[label][c.name]=mu
                rows.append(dict(label=label,model=c.name,**locked(x[good],rr[good]),deviance=deviance(y[good],mu[good]),residual_rms=np.sqrt(np.mean(rr[good]**2)),lag1=np.corrcoef(rr[:-1][good[:-1]&good[1:]],rr[1:][good[:-1]&good[1:]])[0,1],edf='adaptive/unreported',fit_bins=r['background_fit_bins']))
            except Exception as e:rows.append(dict(label=label,model=c.name,error=str(e)))
        for name in names:
            try:
                f=fit_model(name,m,y,width);mu=f.mean;means[label][name]=mu;rr=(y-mu)/np.sqrt(mu)
                rows.append(dict(label=label,model=name,**locked(x[good],rr[good]),deviance=deviance(y[good],mu[good]),residual_rms=np.sqrt(np.mean(rr[good]**2)),edf=f.edf,converged=f.converged))
                if name.startswith('pois_poly') or name.startswith('pspline_16_'):
                    profiles.append(dict(label=label,model=name,**profile_locked(name,m,y,width)))
            except Exception as e:rows.append(dict(label=label,model=name,error=str(e)))
        print('fits',label,flush=True)
    write(out/'backgrounds.csv',rows);write(out/'profiles.csv',profiles)
    np.savez_compressed(out/'background_means.npz',**{l+'__'+n:v for l,d in means.items() for n,v in d.items()})
    # Alternative blocked geometries. This is a sensitivity grid, not a new primary selector.
    cvs=[]
    for folds in (3,5,7):
        for block in (1,2,4,8,16,24,32):
            cfg=dataclasses.replace(config,n_folds=folds,block_size=block)
            rr=pair_cv_scores(m,specs['H2']['count'],m,specs['G1']['count'],frozen_background_candidates(),cfg)
            best=min(rr,key=lambda r:r['combined_total_poisson_deviance'])['name']
            cvs.extend(dict(folds=folds,block=block,winner=best,**r) for r in rr)
    write(out/'cv_geometry.csv',cvs);print('geometry complete',flush=True)
    # Counts are thinned, not mass blocks withheld: estimate prediction at observed coordinates.
    # Five overlapping randomized splits are NOT five independent experiments.
    rng=np.random.default_rng(20260906);cv=[]
    for label in ('H2','G1'):
        y=specs[label]['count']
        for split in range(5):
            train=rng.binomial(y.astype(int),.5).astype(float);val=y-train
            for name in names+[c.name for c in historical]:
                try:
                    if name in names:mu=fit_model(name,m,train,width).mean
                    else:mu=evaluate_candidate(m,train,next(c for c in historical if c.name==name),config)[1]
                    cv.append(dict(label=label,split=split,model=name,predictive_deviance=deviance(val[good],mu[good])))
                except Exception as e:cv.append(dict(label=label,split=split,model=name,error=str(e)))
    write(out/'thinned_predictive.csv',cv)
    # Independent file transfer; one normalization profiled on target eligible counts.
    transfer=[]
    for source,target in [('H1','H2'),('H2','G1'),('G1','G2')]:
        for name,mu in means[source].items():
            y=specs[target]['count'];scale=y[good].sum()/mu[good].sum()
            transfer.append(dict(source=source,target=target,model=name,scale=scale,predictive_deviance=deviance(y[good],scale*mu[good]),**locked(x[good],((y-scale*mu)/np.sqrt(scale*mu))[good])))
    write(out/'independent_transfer.csv',transfer)
    # Rebin selected masses, including shifts, quantile bins and different linear/log grids.
    robust=[];frequency=[];influence=[]
    for label in ('H2','G1','G2'):
        masses=np.load(base/(label+'_masses.npz'))['masses']
        variants=[]
        for n in (175,250,300,350,400,500,700):
            for shift in (0.,.25,.5,.75):
                edges=np.exp(np.linspace(np.log(2),np.log(120),n+1)+shift*np.log(60)/n)
                variants.append((f'log_{n}_shift{shift}',edges))
        for n in (175,350,700):variants.append((f'linear_{n}',np.linspace(2,120,n+1)))
        for n in (175,350,500):variants.append((f'quantile_{n}',np.quantile(masses[(masses>=2)&(masses<=120)],np.linspace(0,1,n+1))))
        for lo,hi in ((2.2,110),(2.5,110),(3,120),(2,70),(4,70),(12,70),(2,30)):
            variants.append((f'window_{lo}_{hi}',np.geomspace(lo,hi,351)))
        for name,edges in variants:
            z=np.sqrt(edges[1:]*edges[:-1]);y=np.histogram(masses,edges)[0].astype(float)
            for model in ('cheb_d7','cheb_d12','spline_s1','spline_s2'):
                try:
                    r,_=evaluate_candidate(z,y,next(c for c in historical if c.name==model),config)
                    robust.append(dict(label=label,variant=name,model=model,**r))
                except Exception as e:robust.append(dict(label=label,variant=name,model=model,error=str(e)))
        y=specs[label]['count']
        # Actual deletion from the background training AND spectral test, not zeroing a bin.
        for size in (1,8,16,32):
            starts=np.flatnonzero(good) if size==1 else range(0,len(m)-size+1,8)
            for start in starts:
                keep=np.ones(len(m),bool);keep[start:start+size]=False
                try:
                    r,_=evaluate_candidate(m[keep],y[keep],BackgroundCandidate('cheb_d7','chebyshev',7),config)
                    influence.append(dict(label=label,start=int(start),size=size,**r))
                except Exception as e:influence.append(dict(label=label,start=int(start),size=size,error=str(e)))
        for name in ('cheb_d7','cheb_d12','spline_s1','spline_s2','pois_poly12'):
            mu=means[label][name];rr=((y-mu)/np.sqrt(mu))[good]
            for low,hi,steps in ((.5,80,3000),(3.1,80,1000),(5,9,1001),(.5,160,4000)):
                _,best=scan_omegas(x[good],rr,np.linspace(low,hi,steps))
                frequency.append(dict(label=label,model=name,search=f'{low}:{hi}:{steps}',**dataclasses.asdict(best)))
            for w in (OMEGA/2,OMEGA,OMEGA*2,OMEGA-1,OMEGA+1,9.7):
                r=weighted_linear_sinusoid(x[good],rr,w)
                frequency.append(dict(label=label,model=name,search='fixed_control',**dataclasses.asdict(r)))
            # Periodic in mass: a distinct template class with 2--50 cycles across domain.
            _,best=scan_omegas(m[good],rr,np.linspace(2,50,1000)*2*np.pi/118)
            frequency.append(dict(label=label,model=name,search='linear_mass_2to50cycles',**dataclasses.asdict(best)))
    write(out/'binning_windows.csv',robust);write(out/'deletions.csv',influence);write(out/'frequency.csv',frequency)
    print('complete',flush=True)

if __name__=='__main__':main()
