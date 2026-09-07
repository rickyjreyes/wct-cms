"""Retrospective predictive stacking; weights never use waveform scores."""
import pathlib,json
import numpy as np
from scipy.optimize import minimize
from models import *

def main():
    base=pathlib.Path('audit/2026-09-06/baseline');out=base.parent/'attacks';names=candidates();rng=np.random.default_rng(20260906);pred=[];targets=[]
    for label in ('H2','G1'):
        s=np.genfromtxt(base/label/'spectrum.csv',delimiter=',',names=True);m=s['mass_center'];width=s['mass_high']-s['mass_low'];g=eligible(m);y=s['count']
        for _ in range(5):
            train=rng.binomial(y.astype(int),.5).astype(float);pred.append(np.column_stack([fit_model(name,m,train,width).mean[g] for name in names]));targets.append((y-train)[g])
    M=np.vstack(pred);y=np.concatenate(targets)
    def fun(w):
        mu=M@w;return np.sum(mu-y*np.log(mu)),M.T@(1-y/mu)
    fit=minimize(fun,np.ones(len(names))/len(names),jac=True,bounds=[(0,1)]*len(names),constraints={'type':'eq','fun':lambda w:w.sum()-1,'jac':lambda w:np.ones(len(w))},method='SLSQP',options={'ftol':1e-8,'maxiter':500})
    if not fit.success:raise RuntimeError(fit.message)
    means=np.load(out/'background_means.npz');evaluation={}
    for label in ('H2','G1','G2'):
        y=np.genfromtxt(base/label/'spectrum.csv',delimiter=',',names=True)['count'];mu=sum(v*means[label+'__'+n] for n,v in zip(names,fit.x));evaluation[label]=dict(**locked(np.log(m[g]),((y-mu)/np.sqrt(mu))[g]),deviance=deviance(y[g],mu[g]))
    (out/'model_averaging.json').write_text(json.dumps(dict(weights=dict(zip(names,map(float,fit.x))),evaluation=evaluation,scope='retrospective predictive weights; not Bayesian probabilities or calibrated significance'),indent=2))
if __name__=='__main__':main()
