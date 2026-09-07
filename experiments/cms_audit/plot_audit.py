"""Render archived numerical evidence; no generated imagery or refitting."""
import csv,json,pathlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=pathlib.Path('audit/2026-09-06')
rows=list(csv.DictReader(open(root/'attacks/backgrounds.csv')))
fig,axes=plt.subplots(1,3,figsize=(14,4.5),layout='constrained')
models=['cheb_d7','cheb_d12','spline_s2','spline_s1','pois_poly12']
for label in ['H2','G1','G2']:
    values=[float(next(r for r in rows if r['label']==label and r['model']==m)['q']) for m in models]
    axes[0].plot(range(len(models)),values,'o-',label=label)
axes[0].set(xticks=range(len(models)),xticklabels=['Cheb 7','Cheb 12','Spline 2','Spline 1','Poisson 12'],ylabel='Positive locked score',title='Background choice dominates')
axes[0].tick_params(axis='x',rotation=30);axes[0].legend()
for case,label in [('historical_null','Historical'),('spline_s1','Spline 1'),('pois_poly12','Poisson 12'),('pspline_16_0.1','Penalized spline')]:
    values=np.array([r['score'] for p in sorted((root/'certificate').glob(case+'_*.json')) for r in json.load(open(p))]);axes[1].hist(values,bins=np.linspace(0,180,55),density=True,histtype='step',label=label)
axes[1].axvline(108.43746213266161,color='black',ls='--',label='Observed')
axes[1].set(title='Same pipeline, different smooth nulls',xlabel='Pair score',ylabel='Probability density',ylim=(0,.045));axes[1].legend(fontsize=8)
axes[1].text(.03,.96,'Historical spike at zero extends above axis',transform=axes[1].transAxes,va='top',fontsize=8)
recovery=list(csv.DictReader(open(root/'controls/recovery_summary.csv')))
for degree in ['7','12']:
    selected=[next(r for r in recovery if r['degree']==degree and r['name']==name) for name in ['amplitude_0.25','amplitude_0.5','nominal']]
    axes[2].errorbar([float(r['amplitude']) for r in selected],[float(r['mean_amplitude']) for r in selected],yerr=[float(r['sd']) for r in selected],fmt='o-',capsize=3,label='Degree '+degree)
axes[2].plot([0,1],[0,1],':',color='black');axes[2].set(title='Flexible fits also absorb real injections',xlabel='Injected amplitude',ylabel='Mean recovered amplitude ± SD');axes[2].legend()
fig.savefig(root/'audit_summary.png',dpi=160)
fig.savefig(root/'audit_summary.pdf')
