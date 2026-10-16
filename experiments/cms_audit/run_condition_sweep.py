"""Frozen retrospective detector-condition sweep; never a discovery search."""
from pathlib import Path
import argparse, csv, hashlib, json, platform, sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'src'),str(ROOT/'scripts'),str(Path(__file__).parent)]
from trace_cms_verified_root_stages import published_identity,hash_input
from cms_wct.cmsio import load_golden_json
from cms_wct.background_cv import BackgroundCandidate,blocked_fold_ids,fit_candidate_train_mask
from cms_wct.background_kill import PairKillConfig,evaluate_candidate,injected_mean
from models import OMEGA,PHASE,MASKS,eligible,fit_model,deviance,locked
from event_samples import extract,selections

MODELS=['cheb_d7','cheb_d12','pspline_cv']
SPLINES=['pspline_16_0.1','pspline_16_1','pspline_16_10','pspline_16_100']
PROTOCOL={
 'schema':'CMS_CONDITION_SWEEP_V1','omega':OMEGA,'phase':PHASE,
 'classification':'retrospective controls, no independent detector or discovery claim',
 'labels':['H1','H2','G1','G2'],'samples':['early','late','random'],
 'entry_cap_before_certification':100000,'seed':20261015,
 'selections':['baseline32','float64','same_sign','one_pair_events','pt5','pt6',
  'eta2p1','medium','isolation0p15','small_dxy','all supported historical muon HLT bits',
  'entry_parity0','entry_parity1','barrel_eta_le_0.9','forward_eta_ge_1.5',
  'min_pt_4_to_6','pt_plateau_proxy_20_10','each run',
  'each run split by floor(lumi/50) modulo 2'],
 'mass_variants':{'nominal':[2,120,350],'lower5':[5,120,350],'expanded_masks':[2,120,350]},
 'mask_windows':MASKS,'expanded_mask_width_factor':1.15,
 'mass_variant_scope':'nominal for every selection; lower5 and expanded_masks for float64 only',
 'backgrounds':MODELS,'spline_candidates':SPLINES,
 'predictive_cv':{'folds':5,'block_bins':16,'selection':'minimum held-out Poisson deviance, no target score'},
 'low_count_rule':'retain row; do not fit if in-range pairs < 1000',
 'uncertainty':'event-cluster Poisson sandwich SE, conditional on fitted background and analysis bins; no background/refit uncertainty',
 'injection':'deterministic +1 Pearson-unit waveform; refit each fixed background, subtract uninjected refit; not power',
 'systematic_experiment':{'trials':200,'seed':20261016,
   'generator':'shared pooled early H2/G1/G2 smooth Poisson spline, CV-selected without WCT scores; scale by each file count',
   'statistic':'minimum positive-only locked Q over H2/G1/G2 for degree7 and degree12',
   'scope':'conditional fitted-smooth-shape simulation; not an instrument-response model or independent null calibration'},
 'interpretation':{
  'stable_positive':'all evaluable file/sample baseline amplitudes positive under a given background',
  'selection_dependent':'any sign reversal among evaluable partitions; report explicitly with counts and SE',
  'model_dependent':'flexible background reduces baseline Q or reverses amplitude; report all, do not select significance',
  'low_power':'in-range <1000, fewer than 8 analysis bins, or conditional |amplitude|/SE < 2',
  'multiplicity':'overlapping controls are descriptive; no minimum p or post-hoc discovery threshold',
  'causality':'localization cannot identify detector cause without efficiency/trigger-object calibration'}
}

def protocol_bytes():
 return (json.dumps(PROTOCOL,indent=2,sort_keys=True)+'\n').encode()

def control_masks(d):
 out=selections(d);b=d['os']&d['tight']
 out.update(barrel_eta_le_0p9=(b&(d['max_eta']<=.9),'mass64'),
  forward_eta_ge_1p5=(b&(d['max_eta']>=1.5),'mass64'),
  min_pt_4_to_6=(b&(d['min_pt']<6),'mass64'),
  pt_plateau_proxy_20_10=(b&(d['max_pt']>=20)&(d['min_pt']>=10),'mass64'))
 for run in sorted(np.unique(d['run'])):
  r=b&(d['run']==run);out[f'run_{run}']=(r,'mass64')
  for parity in (0,1):out[f'run_{run}_lumi50_parity{parity}']=(r&((d['lumi']//50)%2==parity),'mass64')
 return out

def fit(name,m,y,w,good,train=None):
 if name.startswith('cheb'):
  candidate=BackgroundCandidate(name,'chebyshev',int(name.split('d')[-1]))
  if train is None:
   _,mu=evaluate_candidate(m,y,candidate,PairKillConfig(OMEGA,PHASE,tuple())) if good.all() else (None,None)
   if mu is not None:return mu
   # Same robust fitter used by the historical runner, with exact mask intervals.
   from cms_wct.background import robust_smooth_background
   # Explicit training mask uses the identical clipped log-count algorithm.
   return fit_candidate_train_mask(m,y,candidate,good)
  return fit_candidate_train_mask(m,y,candidate,train)
 f=fit_model(name,m,y,w,good if train is None else train)
 if not f.converged:raise RuntimeError('Poisson spline did not converge')
 return f.mean

def cv_score(name,m,y,w,good):
 folds=blocked_fold_ids(good,n_folds=5,block_size=16);total=0.;n=0
 for fold in range(5):
  test=folds==fold;train=good&~test
  if not test.any():continue
  mu=fit(name,m,y,w,good,train)
  total+=deviance(y[test],mu[test]);n+=int(test.sum())
 if n==0:raise RuntimeError('no predictive validation bins')
 return total/n

def statistic(m,y,mu,good):
 mask=good&(mu>=5)
 if mask.sum()<8:raise RuntimeError('fewer than 8 analysis bins')
 r=(y[mask]-mu[mask])/np.sqrt(np.maximum(mu[mask],1))
 return locked(np.log(m[mask]),r),mask

def cluster_se(m,mu,mask,bin_idx,event_idx):
 g=np.cos(OMEGA*np.log(m[mask])-PHASE);g-=g.mean()
 influence=np.zeros(len(m));influence[mask]=g/(g@g)/np.sqrt(np.maximum(mu[mask],1))
 _,ix=np.unique(event_idx,return_inverse=True)
 scores=np.bincount(ix,weights=influence[bin_idx])
 return float(np.sqrt(scores@scores))

def evaluate(d,selected,mass_name,variant,label,sampling,control):
 lo,hi,bins=PROTOCOL['mass_variants'][variant];edges=np.geomspace(lo,hi,bins+1)
 m=np.sqrt(edges[:-1]*edges[1:]);w=np.diff(edges)
 mass=d[mass_name][selected];events=d['entry'][selected]
 bi=np.searchsorted(edges,mass,side='right')-1;bi[mass==hi]=bins-1
 valid=(bi>=0)&(bi<bins)&np.isfinite(mass);bi=bi[valid];events=events[valid]
 y=np.bincount(bi,minlength=bins).astype(float)
 windows=MASKS if variant!='expanded_masks' else tuple(((a+b)/2-.575*(b-a),(a+b)/2+.575*(b-a)) for a,b in MASKS)
 good=eligible(m,windows)
 common=dict(label=label,sampling=sampling,control=control,mass_variant=variant,
  selected_pairs=int(selected.sum()),in_range_pairs=len(bi),event_count=len(np.unique(events)),
  events_unit='source ROOT entry after certification; raw-key uniqueness audited',
  uncertainty_method='event-cluster Poisson SE conditional on background',
  trigger_scope='HLT bit / kinematic proxy; no trigger-object matching or efficiency unfolding')
 cv={};errors={}
 if len(bi)>=1000:
  for name in ['cheb_d7','cheb_d12']+SPLINES:
   try:cv[name]=cv_score(name,m,y,w,good)
   except (ValueError,RuntimeError,np.linalg.LinAlgError) as exc:errors[name]=str(exc)
 possible=[name for name in SPLINES if name in cv and np.isfinite(cv[name])]
 winner=min(possible,key=lambda n:cv[n]) if possible else None
 rows=[]
 for model in MODELS:
  name=winner if model=='pspline_cv' else model
  row={**common,'model':model,'fitted_model':name,'predictive_deviance_per_bin':cv.get(name),
       'spline_cv_candidates':json.dumps({k:cv.get(k,errors.get(k)) for k in SPLINES}) if model=='pspline_cv' else ''}
  if len(bi)<1000:row.update(status='LOW_COUNT',low_power=True)
  elif name is None:row.update(status='FIT_FAILED',error=json.dumps(errors),low_power=True)
  else:
   try:
    mu=fit(name,m,y,w,good);score,mask=statistic(m,y,mu,good)
    se=cluster_se(m,mu,mask,bi,events)
    mu0=fit(name,m,mu,w,good);s0,_=statistic(m,mu,mu0,good)
    injected=mu.copy();injected[mask]+=np.sqrt(mu[mask])*np.cos(OMEGA*np.log(m[mask])-PHASE)
    injected=np.maximum(injected,1e-9)
    mui=fit(name,m,injected,w,good);si,_=statistic(m,injected,mui,good)
    row.update(status='EVALUATED',amplitude=score['amplitude'],positive_sign=score['amplitude']>0,
     Q=score['q'],analysis_bins=int(mask.sum()),conditional_cluster_se=se,
     conditional_ci95_lower=score['amplitude']-1.96*se,conditional_ci95_upper=score['amplitude']+1.96*se,
     low_power=se<=0 or abs(score['amplitude'])<2*se,
     deterministic_injection_retention=si['amplitude']-s0['amplitude'],
     injection_power=None,fit_deviance_per_bin=deviance(y[good],mu[good])/good.sum())
   except (ValueError,RuntimeError,np.linalg.LinAlgError) as exc:row.update(status='FIT_FAILED',error=str(exc),low_power=True)
  rows.append(row)
 return rows,(m,y,w,good)

def systematic_experiment(spectra):
 m=spectra['H2'][0];w=spectra['H2'][2];good=spectra['H2'][3]
 pooled=sum(spectra[k][1] for k in ['H2','G1','G2']);scores={}
 for name in SPLINES:
  try:scores[name]=cv_score(name,m,pooled,w,good)
  except (ValueError,RuntimeError,np.linalg.LinAlgError):pass
 if not scores:return {'status':'FIT_FAILED','interpretation':'shared-shape experiment not executed'}
 best=min(scores,key=scores.get);shape=fit(best,m,pooled,w,good)
 rng=np.random.default_rng(PROTOCOL['systematic_experiment']['seed']);obs={};trials=[]
 for name in ['cheb_d7','cheb_d12']:
  obs[name]=min(statistic(m,spectra[k][1],fit(name,m,spectra[k][1],w,good),good)[0]['q'] for k in ['H2','G1','G2'])
 for trial in range(PROTOCOL['systematic_experiment']['trials']):
  pseudo={k:rng.poisson(shape*(spectra[k][1].sum()/pooled.sum())) for k in ['H2','G1','G2']}
  for name in ['cheb_d7','cheb_d12']:
   try:
    a=[statistic(m,y,fit(name,m,y,w,good),good)[0] for y in pseudo.values()]
    trials.append(dict(trial=trial,model=name,min_Q=min(s['q'] for s in a),all_positive=all(s['amplitude']>0 for s in a)))
   except (ValueError,RuntimeError,np.linalg.LinAlgError) as exc:trials.append(dict(trial=trial,model=name,error=str(exc)))
 return dict(status='EXECUTED',generator=best,predictive_cv=scores,observed_min_Q=obs,trials=trials,
  scope=PROTOCOL['systematic_experiment']['scope'],calibrated_discovery_p=None)

def main():
 p=argparse.ArgumentParser();p.add_argument('--freeze',type=Path);p.add_argument('--protocol',type=Path)
 p.add_argument('--input-root',type=Path);p.add_argument('--output',type=Path);a=p.parse_args()
 if a.freeze:
  a.freeze.parent.mkdir(parents=True,exist_ok=True)
  with a.freeze.open('xb') as f:f.write(protocol_bytes())
  print(hashlib.sha256(protocol_bytes()).hexdigest());return
 if not a.protocol or a.protocol.read_bytes()!=protocol_bytes():p.error('protocol bytes differ from frozen implementation')
 if not a.input_root or not a.output:p.error('input root and new output directory required')
 a.output.mkdir(parents=True,exist_ok=False)
 cert=load_golden_json(str(ROOT/'data/cert/14221/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON_MuonPhys.txt'))
 rows=[];cutflows=[];receipts=[];spectra={};raw_keys={}
 for label in PROTOCOL['labels']:
  expected=published_identity(ROOT,label);path=a.input_root/expected['filename']
  receipts.append(dict(label=label,**hash_input(path,expected)))
  for sampling in PROTOCOL['samples']:
   d,flow,keys=extract(path,cert,sampling,PROTOCOL['entry_cap_before_certification'],PROTOCOL['seed'])
   cutflows.append(dict(label=label,**flow));raw_keys[label+'_'+sampling]=keys
   for name,(selected,mass) in control_masks(d).items():
    variants=['nominal','lower5','expanded_masks'] if name=='float64' else ['nominal']
    for variant in variants:
     out,spectrum=evaluate(d,selected,mass,variant,label,sampling,name);rows.extend(out)
     if sampling=='early' and name=='float64' and variant=='nominal':spectra[label]=spectrum
   print(label,sampling,len(rows),flush=True)
   (a.output/'rows_checkpoint.json').write_text(json.dumps(rows,indent=2)+'\n')
 keys=sorted({k for row in rows for k in row})
 with (a.output/'results.csv').open('x') as f:
  writer=csv.DictWriter(f,fieldnames=keys);writer.writeheader();writer.writerows(rows)
 systematic=systematic_experiment(spectra)
 (a.output/'common_shape_simulation.json').write_text(json.dumps(systematic,indent=2)+'\n')
 report=dict(schema='CMS_CONDITION_SWEEP_EXECUTION_V1',protocol_sha256=hashlib.sha256(a.protocol.read_bytes()).hexdigest(),
  source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),Path(__file__).with_name('event_samples.py'),Path(__file__).with_name('models.py')]},
  environment=dict(python=platform.python_version(),numpy=np.__version__),input_receipts=receipts,cutflows=cutflows,
  rows=len(rows),overlaps={a+'__'+b:len(raw_keys[a]&raw_keys[b]) for a in raw_keys for b in raw_keys if a<b},
  scientific_classification='RETROSPECTIVE_DESCRIPTIVE_PENDING_REVIEW',independent_detector_replication=False,
  limitations=['Conditional uncertainty excludes background refit/selection','Injection retention is not detection power','HLT bits and pT proxy do not calibrate trigger efficiency','Shared-shape simulation does not identify detector causality'])
 (a.output/'execution.json').write_text(json.dumps(report,indent=2)+'\n')
 print('completed',len(rows),flush=True)
if __name__=='__main__':main()
