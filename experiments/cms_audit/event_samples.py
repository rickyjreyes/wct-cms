"""Frozen entry sampling and event/pair provenance shared by the audit runners."""
from pathlib import Path
import hashlib
import awkward as ak
import numpy as np
import uproot
from cms_wct.cmsio import golden_mask, invariant_mass_from_pairs
FIELDS=['pt','eta','phi','mass','charge','tightId','mediumId','pfRelIso04_all','dxy']

def select_entries(n, sampling, cap=100000, seed=20261015):
    count=min(n,cap)
    if sampling=='early': return np.arange(count,dtype=np.int64)
    if sampling=='late': return np.arange(n-count,n,dtype=np.int64)
    if sampling=='random':return np.sort(np.random.default_rng(seed).choice(n,count,replace=False))
    raise ValueError('unknown sampling policy')

def key_digest(keys):
    return hashlib.sha256(np.asarray(sorted(keys),dtype='<u8').reshape(-1,3).tobytes()).hexdigest()

def extract(path, cert, sampling='early', cap=100000, seed=20261015):
    with uproot.open(path) as f:
        t=f['Events'];entries=select_entries(t.num_entries,sampling,cap,seed)
        triggers=[k for k in t.keys() if k.startswith('HLT_') and
                  ('Mu17_TrkIsoVVL_Mu8_TrkIsoVVL' in k or k in ['HLT_DoubleMu0','HLT_Mu17_Mu8_DZ'])]
        fields=['run','luminosityBlock','event']+['Muon_'+k for k in FIELDS]+triggers
        chunks=[]
        for block in np.unique(entries//100000):
            selected=entries[entries//100000==block]
            start=int(selected[0]);stop=int(selected[-1])+1
            a=t.arrays(fields,entry_start=start,entry_stop=stop,how=dict)
            chunks.append({k:v[selected-start] for k,v in a.items()})
        if not chunks:raise ValueError('empty ROOT input')
        ar={k:ak.concatenate([c[k] for c in chunks]) for k in fields}
        total=int(t.num_entries)
    run=ak.to_numpy(ar['run']);lumi=ak.to_numpy(ar['luminosityBlock']);ev=ak.to_numpy(ar['event'])
    rawkeys=set(zip(run.tolist(),lumi.tolist(),ev.tolist()))
    gold=golden_mask(run,lumi,cert);ar={k:v[gold] for k,v in ar.items()}
    certkeys=set(zip(run[gold].tolist(),lumi[gold].tolist(),ev[gold].tolist()))
    keep=(ar['Muon_pt']>=4)&(abs(ar['Muon_eta'])<=2.4)
    mu=ak.zip({k:ar['Muon_'+k][keep] for k in FIELDS})
    pairs=ak.combinations(mu,2,fields=['a','b'])
    flat=lambda a:ak.to_numpy(ak.flatten(a))
    broadcast=lambda a:flat(ak.broadcast_arrays(a,pairs.a.pt)[0])
    out={'entry':broadcast(ak.local_index(ar['run'])).astype(np.int64),
         'source_entry':broadcast(entries[gold]).astype(np.int64),
         'run':broadcast(ar['run']), 'lumi':broadcast(ar['luminosityBlock']), 'event':broadcast(ar['event']),
         'mass32':flat(invariant_mass_from_pairs(pairs))}
    ma=ak.zip({k:ak.values_astype(pairs.a[k],np.float64) for k in ['pt','eta','phi','mass']})
    mb=ak.zip({k:ak.values_astype(pairs.b[k],np.float64) for k in ['pt','eta','phi','mass']})
    out['mass64']=flat(invariant_mass_from_pairs(ak.zip({'a':ma,'b':mb})))
    out.update(os=flat(pairs.a.charge*pairs.b.charge<0),tight=flat(pairs.a.tightId & pairs.b.tightId),
       medium=flat(pairs.a.mediumId & pairs.b.mediumId),
       min_pt=flat(np.minimum(pairs.a.pt,pairs.b.pt)),max_pt=flat(np.maximum(pairs.a.pt,pairs.b.pt)),
       max_eta=flat(np.maximum(abs(pairs.a.eta),abs(pairs.b.eta))),
       max_iso=flat(np.maximum(pairs.a.pfRelIso04_all,pairs.b.pfRelIso04_all)),
       max_dxy=flat(np.maximum(abs(pairs.a.dxy),abs(pairs.b.dxy))))
    for tr in triggers:out[tr]=broadcast(ar[tr]).astype(bool)
    cutflow={'root_entries':total,'sampled_before_certification':len(entries),'certified':int(gold.sum()),
      'certification_rejected':int((~gold).sum()),'raw_unique_keys':len(rawkeys),'certified_unique_keys':len(certkeys),
      'raw_event_keys_sha256':key_digest(rawkeys),'certified_event_keys_sha256':key_digest(certkeys),
      'sampled_entry_indices_sha256':hashlib.sha256(entries.astype('<i8').tobytes()).hexdigest(),
      'candidate_pairs':len(out['entry']),'baseline_pairs':int((out['os']&out['tight']).sum()),
      'sampling':sampling,'seed':seed,'cap':cap,
      'trigger_object_matching':'NOT_PERFORMED; HLT bits only, no efficiency unfolding',
      'runs':{str(k):int(v) for k,v in zip(*np.unique(run[gold],return_counts=True))}}
    return out,cutflow,rawkeys

def selections(d):
    base=d['os']&d['tight']; ev,counts=np.unique(d['entry'][base],return_counts=True)
    out={'baseline32':(base,'mass32'),'float64':(base,'mass64'),
      'same_sign':(~d['os']&d['tight'],'mass64'),
      'one_pair_events':(base&np.isin(d['entry'],ev[counts==1]),'mass64'),
      'pt5':(base&(d['min_pt']>=5),'mass64'),'pt6':(base&(d['min_pt']>=6),'mass64'),
      'eta2p1':(base&(d['max_eta']<=2.1),'mass64'),'medium':(d['os']&d['medium'],'mass64'),
      'isolation0p15':(base&(d['max_iso']<.15),'mass64'),'small_dxy':(base&(d['max_dxy']<.02),'mass64')}
    for k in d:
        if k.startswith('HLT_'):out[k]=(base&d[k],'mass64')
    for parity in (0,1):out[f'entry_parity{parity}']=(base&(d['entry']%2==parity),'mass64')
    return out
