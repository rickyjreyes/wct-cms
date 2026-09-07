"""Historical source reproduction from fetch_inputs.py downloads; separate output directory."""
import json,pathlib,sys,dataclasses,numpy as np
from cms_wct.cli import parser
from cms_wct.analysis import run_analysis
from cms_wct.cmsio import extract_dimuon_masses,load_golden_json
from cms_wct.background_kill import evaluate_pair_pipeline,PairKillConfig
sys.path.insert(0,'scripts')
import run_phase_locked_period as locked
root=pathlib.Path('data/audit_inputs')
out=pathlib.Path('results/cms_audit_clean');out.mkdir(parents=True,exist_ok=True)
cert='data/cert/14221/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON_MuonPhys.txt'
base=['--golden-json',cert,'--mass-min','2','--mass-max','120','--bins','350','--log-bins','--tight-id','--max-events','100000','--step-size','100 MB','--fit-degree','7']
for w in ['2.9:3.3','3.55:3.85','8.5:11.5','80:100']:base+=['--mask-window',w]
for label,prefix in [('H1','127C'),('H2','183B'),('G1','05DD'),('G2','209D')]:
 path=str(next(root.glob(prefix+'*.root')))
 args=parser().parse_args(base+['--input',path,'--output-dir',str(out/label),'--omega-min','3.1','--omega-max','80','--omega-steps','1000','--frozen-omega','7.025825825825827','--permutations','1000','--parametric-bootstrap','500'])
 masses,counters=extract_dimuon_masses([path],args,load_golden_json(cert))
 np.savez_compressed(out/(label+'_masses.npz'),masses=masses)
 result=run_analysis(args)
 print(label,json.dumps(dataclasses.asdict(result)),flush=True)
 if label=='G2':
  sys.argv=['run_phase_locked_period.py']+base+['--input',path,'--output-dir',str(out/'G2_locked'),'--omega','7.025825825825827','--phase','-0.2313916852932179']
  locked.main()
h=np.genfromtxt(out/'H2/spectrum.csv',delimiter=',',names=True)
g=np.genfromtxt(out/'G1/spectrum.csv',delimiter=',',names=True)
c=PairKillConfig(omega=7.025825825825827,phase=-.2313916852932179,excluded_windows=((2.9,3.3),(3.55,3.85),(8.5,11.5),(80,100)))
r,*_=evaluate_pair_pipeline(h['mass_center'],h['count'],g['mass_center'],g['count'],None,c)
(out/'pair_observed.json').write_text(json.dumps(r,indent=2))
print('PAIR',r['selected_background'],r['primary_pair_score_min_locked_delta_chi2'],flush=True)
