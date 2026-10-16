#!/usr/bin/env python3
"""Preserve prior observations and emit paired evidence, including discrepancies."""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd

def compare(directory, prior_path):
 prior=json.loads(prior_path.read_text());new=json.loads((directory/'summary.json').read_text())
 sel=pd.read_csv(directory/'selections.csv');boot=pd.read_csv(directory/'bootstrap.csv')
 rows=[]
 for old in prior['files']:
  label=old['label'];now=next(f for f in new['files'] if f['label']==label)
  for key in ['certified','pairs','pair_events','multiple_pair_events','max_pairs','runs','float64_count_L1']:
   if key in old:rows.append(dict(label=label,metric=key,prior=old[key],regenerated=now[key],exact_match=old[key]==now[key]))
  for degree in (7,12):
   obs=float(sel.loc[(sel.label==label)&(sel.selection=='same_sign')&(sel.degree==degree),'locked_delta_chi2'].iloc[0])
   expected=old[f'same_sign_degree{degree}_score']
   rows.append(dict(label=label,metric=f'same_sign_degree{degree}_score',prior=expected,regenerated=obs,absolute_delta=obs-expected))
  vals=boot.loc[(boot.label==label)&(boot.degree==7),'amplitude']
  interval=np.quantile(vals,[.025,.975]).tolist()
  rows.append(dict(label=label,metric='degree7_cluster_interval',prior=old['degree7_cluster_interval'],regenerated=interval,trials=len(vals)))
  if 'float64_score' in old:
   obs=float(sel.loc[(sel.label==label)&(sel.selection=='float64')&(sel.degree==7),'locked_delta_chi2'].iloc[0])
   rows.append(dict(label=label,metric='float64_score',prior=old['float64_score'],regenerated=obs,absolute_delta=obs-old['float64_score']))
 return dict(schema='CMS_PRIOR_EVENT_COMPARISON_V1',pairs=rows,raw_overlaps=new['overlaps'],
   counts_match=all(r.get('exact_match',True) for r in rows),
   interpretation='Paired historical reproduction; numerical differences retained, no significance promotion.')

def main():
 p=argparse.ArgumentParser();p.add_argument('--directory',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
 a=p.parse_args();out=compare(a.directory,Path('audit/2026-09-06/prior_pass_event_observations.json'))
 with a.output.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
