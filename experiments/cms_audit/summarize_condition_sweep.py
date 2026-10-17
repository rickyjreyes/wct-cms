"""Report the entire fixed matrix, including failures; no significance selection."""
import argparse,json
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
 p=argparse.ArgumentParser();p.add_argument('--directory',type=Path,required=True);a=p.parse_args()
 d=pd.read_csv(a.directory/'results.csv');b=d[(d.control=='float64')&(d.mass_variant=='nominal')]
 summary={'schema':'CMS_CONDITION_SWEEP_SUMMARY_V1','rows':len(d),'status_counts':d.status.value_counts().to_dict(),
 'classification':'BROAD_DEGREE7_REPEATABILITY_WITH_BACKGROUND_AND_SELECTION_DEPENDENCE',
 'full_sweep_calibrated_p':None,'independent_detector_replication':False,'models':{}}
 for model in ['cheb_d7','cheb_d12','pspline_cv']:
  z=d[(d.model==model)&(d.status=='EVALUATED')];base=b[b.model==model]
  summary['models'][model]={'evaluated_rows':len(z),'nonpositive_amplitudes':int((z.amplitude<=0).sum()),
   'conditional_low_power_rows':int(z.low_power.sum()),'baseline_positive':int((base.amplitude>0).sum()),
   'baseline_amplitude_range':[float(base.amplitude.min()),float(base.amplitude.max())],
   'baseline_Q_range':[float(base.Q.min()),float(base.Q.max())],
   'baseline_injection_retention_range':[float(base.deterministic_injection_retention.min()),float(base.deterministic_injection_retention.max())]}
 sim=json.loads((a.directory/'common_shape_simulation.json').read_text());s={}
 for model,obs in sim['observed_min_Q'].items():
  trials=[t for t in sim['trials'] if t['model']==model];ok=[t for t in trials if 'min_Q' in t]
  s[model]={'observed_min_Q':obs,'evaluated_trials':len(ok),'failed_trials':len(trials)-len(ok),
   'exceedances':sum(t['min_Q']>=obs for t in ok),'all_positive_trials':sum(t['all_positive'] for t in ok),
   'min_Q_quantiles_05_50_95':np.quantile([t['min_Q'] for t in ok],[.05,.5,.95]).tolist()}
 summary['shared_smooth_shape_experiment']=s
 summary['limitations']=['Conditional uncertainty omits background fitting/selection','Deterministic injection retention is not power',
  'Degree12 has poor blocked predictive extrapolation here; spline is the competitive alternative',
  'Common-shape simulation shows repeatability but does not fully reproduce observed degree7 joint score',
  'No trigger-object or response-efficiency calibration; no causal detector identification']
 with (a.directory/'summary.json').open('x') as f:json.dump(summary,f,indent=2);f.write('\n')
 labels=[l+' '+s for l in ['H1','H2','G1','G2'] for s in ['early','late','random']]
 fig,axes=plt.subplots(2,1,figsize=(11,8),sharex=True,layout='constrained')
 colors=['#135c87','#bb4b31','#2f856b']
 for j,(model,color) in enumerate(zip(['cheb_d7','cheb_d12','pspline_cv'],colors)):
  sub=b[b.model==model].copy();sub.index=sub.label+' '+sub.sampling;sub=sub.loc[labels];x=np.arange(12)+(j-1)*.20
  axes[0].errorbar(x,sub.amplitude,yerr=1.96*sub.conditional_cluster_se,fmt='o',ms=4,capsize=2,label=model,color=color)
  axes[1].scatter(x,sub.Q,label=model,color=color,s=25)
 axes[0].axhline(0,color='#999999',lw=.8);axes[0].set_ylabel('Locked amplitude ± 1.96 conditional SE')
 axes[1].set_ylabel('Positive-only locked Q');axes[1].set_xticks(np.arange(12),labels,rotation=45,ha='right')
 axes[0].legend(ncol=3,frameon=False);axes[0].set_title('CMS 2016: broad historical repeatability, strong background dependence',loc='left',pad=14)
 for ax in axes:ax.grid(axis='y',alpha=.2);ax.spines[['top','right']].set_visible(False)
 fig.supxlabel('All 12 fixed file/sample combinations · 100,000 entries before certification\nUncertainty conditional on background; these overlapping same-detector controls are not independent replications.',fontsize=9)
 fig.savefig(a.directory/'condition_summary.svg');fig.savefig(a.directory/'condition_summary.png',dpi=150)
 print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
