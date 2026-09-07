"""Reproduce and stress the historical null, with durable 100-trial checkpoints.

The original selector/refit/statistic are unchanged. Alternative generating
means contain no explicitly injected waveform. Search statistics are explicitly
permissive sensitivity checks, not an exact reconstruction of missing history.
"""
import argparse,concurrent.futures,csv,json,pathlib
import numpy as np
from cms_wct.background_kill import PairKillConfig,evaluate_pair_pipeline,injected_mean
from cms_wct.signature import _prepare_unweighted_basis,_scores_from_basis
from models import OMEGA,PHASE,MASKS,eligible


def evaluate(job):
    case,trial,m,a,b,search=job
    r,ba,bb=evaluate_pair_pipeline(m,a,m,b,None,PairKillConfig(OMEGA,PHASE,MASKS))
    row=dict(case=case,trial=trial,score=r['primary_pair_score_min_locked_delta_chi2'],selected=r['selected_background']['name'],amp_a=r['spectrum_a']['locked_signed_amplitude'],amp_b=r['spectrum_b']['locked_signed_amplitude'])
    if search:
        omegas=np.unique(np.r_[np.linspace(.5,160,2000),OMEGA]);scores=[]
        for y,bg in [(a,ba),(b,bb)]:
            g=eligible(m)&(bg>=5);x=np.log(m[g]);res=((y-bg)/np.sqrt(bg))[g]
            scores.append(_scores_from_basis(res,*_prepare_unweighted_basis(x,omegas)))
        row['broad_score']=float(np.max(np.minimum(*scores)))
    return row


def main():
    p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=4);p.add_argument('--out',default='audit/2026-09-06/certificate');a=p.parse_args();out=pathlib.Path(a.out);out.mkdir(parents=True,exist_ok=True)
    base=pathlib.Path('audit/2026-09-06/baseline');h=np.genfromtxt(base/'H2/spectrum.csv',delimiter=',',names=True);g=np.genfromtxt(base/'G1/spectrum.csv',delimiter=',',names=True);m=h['mass_center'];means=np.load('audit/2026-09-06/attacks/background_means.npz');cfg=PairKillConfig(OMEGA,PHASE,MASKS)
    observed,ma,mb=evaluate_pair_pipeline(m,h['count'],m,g['count'],None,cfg);T=observed['primary_pair_score_min_locked_delta_chi2']
    jobs_by_case={};rng=np.random.default_rng(20260901)
    jobs_by_case['historical_null']=[('historical_null',i,m,rng.poisson(ma).astype(float),rng.poisson(mb).astype(float),i<1000) for i in range(10000)]
    for ci,name in enumerate(['spline_s1','pois_poly12','pspline_16_0.1'],1):
        jobs=[]
        for i in range(1000):
            rng=np.random.default_rng(20260907+ci*100000+i)
            jobs.append((name,i,m,rng.poisson(means['H2__'+name]).astype(float),rng.poisson(means['G1__'+name]).astype(float),False))
        jobs_by_case[name]=jobs
    rng=np.random.default_rng(20260902)
    for amp in (.25,.5,.75,1.):
        name=f'injection_{amp}';aa=injected_mean(m,ma,cfg,amp);bb=injected_mean(m,mb,cfg,amp)
        jobs_by_case[name]=[(name,i,m,rng.poisson(aa).astype(float),rng.poisson(bb).astype(float),False) for i in range(1000)]
    jobs=[]
    for i in range(1000):
        rng=np.random.default_rng(20260909+i);name='search_poly12'
        jobs.append((name,i,m,rng.poisson(means['H2__pois_poly12']).astype(float),rng.poisson(means['G1__pois_poly12']).astype(float),True))
    jobs_by_case[name]=jobs
    (out/'manifest.json').write_text(json.dumps(dict(source_commit='1c580e141eb5be65b81b0ca638da25bc69d2f0cc',observed=T,config=cfg.to_dict(),case_counts={k:len(v) for k,v in jobs_by_case.items()},historical_seed=20260901,injection_seed=20260902,alternative_seed_rule='20260907 + (1 for spline_s1, 2 for pois_poly12, 3 for pspline_16_0.1)*100000 + trial',search_seed_rule='20260909+trial'),indent=2))
    allrows={}
    with concurrent.futures.ProcessPoolExecutor(a.workers) as pool:
        for name,jobs in jobs_by_case.items():
            rows=[]
            for start in range(0,len(jobs),100):
                f=out/f'{name}_{start:05d}.json'
                if not f.exists():
                    chunk=list(pool.map(evaluate,jobs[start:start+100],chunksize=5));tmp=f.with_suffix('.tmp');tmp.write_text(json.dumps(chunk));tmp.replace(f)
                rows.extend(json.load(open(f)));print(name,start,flush=True)
            allrows[name]=rows
    threshold=float(np.quantile([r['score'] for r in allrows['historical_null']],.95));summaries=[]
    for name,rows in allrows.items():
        q=np.array([r['score'] for r in rows]);k=int(np.sum(q>=T));r=dict(case=name,n=len(rows),exceedances=k,add_one_p=(k+1)/(len(rows)+1),q95=float(np.quantile(q,.95)),median_amp_a=float(np.median([r['amp_a'] for r in rows])),median_amp_b=float(np.median([r['amp_b'] for r in rows])))
        scan=[x['broad_score'] for x in rows if 'broad_score' in x]
        if scan:r.update(broad_n=len(scan),broad_exceedances=int(np.sum(np.array(scan)>=T)),broad_add_one_p=(1+int(np.sum(np.array(scan)>=T)))/(1+len(scan)))
        if name.startswith('injection'):
            amp=float(name.split('_')[1]);r.update(median_min_retention=float(np.median([min(x['amp_a'],x['amp_b']) for x in rows])/amp),power=float(np.mean(q>=threshold)),null_threshold=threshold)
        summaries.append(r)
    (out/'summary.json').write_text(json.dumps(summaries,indent=2));print(json.dumps(summaries,indent=2))
if __name__=='__main__':main()
