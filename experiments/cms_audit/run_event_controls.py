"""Regenerate event-key, precision, selection and cluster-bootstrap controls."""
import pathlib,json
import numpy as np
import awkward as ak
import uproot
from cms_wct.cmsio import golden_mask,load_golden_json,invariant_mass_from_pairs
from cms_wct.background import make_histogram
from cms_wct.background_kill import evaluate_candidate,PairKillConfig
from cms_wct.background_cv import BackgroundCandidate
from models import OMEGA,PHASE,MASKS
from run_spectrum_attacks import write


def main():
    root=pathlib.Path('data/audit_inputs');out=pathlib.Path('audit/2026-09-06/events');out.mkdir(exist_ok=True);ids={};summary=[];rows=[];boot=[];rng=np.random.default_rng(20260911)
    cert=load_golden_json('data/cert/14221/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON_MuonPhys.txt')
    for label,prefix in [('H1','127C'),('H2','183B'),('G1','05DD'),('G2','209D')]:
        file=next(root.glob(prefix+'*.root'));t=uproot.open(file)['Events'];fields=['pt','eta','phi','mass','charge','tightId','mediumId','pfRelIso04_all','dxy']
        triggers=[k for k in t.keys() if k.startswith('HLT_') and ('Mu17_TrkIsoVVL_Mu8_TrkIsoVVL' in k or k in ['HLT_DoubleMu0','HLT_Mu17_Mu8_DZ'])]
        ar=t.arrays(['run','luminosityBlock','event']+['Muon_'+f for f in fields]+triggers,entry_stop=100000,how=dict)
        run=ak.to_numpy(ar['run']);lumi=ak.to_numpy(ar['luminosityBlock']);event=ak.to_numpy(ar['event']);ids[label]=set(zip(run.tolist(),lumi.tolist(),event.tolist()));gold=golden_mask(run,lumi,cert);ar={k:v[gold] for k,v in ar.items()}
        keep=(ar['Muon_pt']>=4)&(abs(ar['Muon_eta'])<=2.4);mu=ak.zip({f:ar['Muon_'+f][keep] for f in fields});pairs=ak.combinations(mu,2,fields=['a','b']);flat=lambda v:ak.to_numpy(ak.flatten(v))
        entry=flat(ak.broadcast_arrays(ak.local_index(ar['run']),pairs.a.pt)[0]).astype(int)
        mass32=flat(invariant_mass_from_pairs(pairs));ma=ak.zip({f:ak.values_astype(pairs.a[f],np.float64) for f in ['pt','eta','phi','mass']});mb=ak.zip({f:ak.values_astype(pairs.b[f],np.float64) for f in ['pt','eta','phi','mass']});mass64=flat(invariant_mass_from_pairs(ak.zip({'a':ma,'b':mb})))
        os=flat(pairs.a.charge*pairs.b.charge<0);tight=flat(pairs.a.tightId & pairs.b.tightId);baseline=os&tight;pe,pc=np.unique(entry[baseline],return_counts=True)
        data=dict(mass32=mass32,mass64=mass64,entry=entry,os=os,tight=tight);np.savez_compressed(out/(label+'_pairs.npz'),**data)
        edges=np.geomspace(2,120,351);old=np.load(f'audit/2026-09-06/baseline/{label}_masses.npz')['masses'];h32=np.histogram(mass32[baseline],edges)[0];h64=np.histogram(mass64[baseline],edges)[0]
        summary.append(dict(label=label,certified=int(gold.sum()),unique_events=len(ids[label]),pairs=int(baseline.sum()),pair_events=len(pc),multiple_pair_events=int(np.sum(pc>1)),max_pairs=int(pc.max()),changed_bins_vs_initial=int(np.sum(h32!=np.histogram(old,edges)[0])),float64_count_L1=int(np.abs(h64-h32).sum()),runs={str(k):int(v) for k,v in zip(*np.unique(run[gold],return_counts=True))}))
        sels={'baseline32':(baseline,mass32),'float64':(baseline,mass64),'same_sign':(~os&tight,mass64),'one_pair_events':(baseline&np.isin(entry,pe[pc==1]),mass64),'pt5':(baseline&flat(np.minimum(pairs.a.pt,pairs.b.pt)>=5),mass64),'pt6':(baseline&flat(np.minimum(pairs.a.pt,pairs.b.pt)>=6),mass64),'eta2p1':(baseline&flat(np.maximum(abs(pairs.a.eta),abs(pairs.b.eta))<=2.1),mass64),'medium':(os&flat(pairs.a.mediumId & pairs.b.mediumId),mass64),'isolation0p15':(baseline&flat(np.maximum(pairs.a.pfRelIso04_all,pairs.b.pfRelIso04_all)<.15),mass64),'small_dxy':(baseline&flat(np.maximum(abs(pairs.a.dxy),abs(pairs.b.dxy))<.02),mass64)}
        for tr in triggers:sels[tr]=(baseline&flat(ak.broadcast_arrays(ar[tr],pairs.a.pt)[0]),mass64)
        for parity in (0,1):sels[f'entry_parity{parity}']=(baseline&(entry%2==parity),mass64)
        for name,(selected,mass) in sels.items():
            y,_,m=make_histogram(mass[selected],2,120,350,True)
            for d in (7,12):
                row=dict(label=label,selection=name,degree=d,in_range=int(y.sum()))
                if y.sum()<1000:row['not_evaluated']='low statistics'
                else:
                    try:r,_=evaluate_candidate(m,y,BackgroundCandidate(f'cheb_d{d}','chebyshev',d),PairKillConfig(OMEGA,PHASE,MASKS));row.update(r)
                    except Exception as e:row['error']=str(e)
                rows.append(row)
        if label!='H1':
            mass=mass32[baseline];ev=entry[baseline];bi=np.searchsorted(edges,mass,side='right')-1;valid=(bi>=0)&(bi<350);bi=bi[valid];ev=ev[valid];m=np.sqrt(edges[1:]*edges[:-1])
            for trial in range(300):
                weights=rng.poisson(1.,int(ev.max())+1);y=np.bincount(bi,weights=weights[ev],minlength=350)
                for d in (7,12):
                    r,_=evaluate_candidate(m,y,BackgroundCandidate(f'cheb_d{d}','chebyshev',d),PairKillConfig(OMEGA,PHASE,MASKS));boot.append(dict(label=label,trial=trial,degree=d,amplitude=r['locked_signed_amplitude']))
        print(label,flush=True)
    (out/'summary.json').write_text(json.dumps(dict(files=summary,overlaps={a+'__'+b:len(ids[a]&ids[b]) for a in ids for b in ids if a<b}),indent=2));write(out/'selections.csv',rows);write(out/'bootstrap.csv',boot)
if __name__=='__main__':main()
