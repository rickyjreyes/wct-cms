import sys
from pathlib import Path
import numpy as np
import pytest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments/cms_audit'))
import event_samples as samples
import run_condition_sweep as sweep

def test_entry_policies_are_frozen_and_pre_certification():
 assert samples.select_entries(150000,'early')[-1]==99999
 assert samples.select_entries(150000,'late')[0]==50000
 a=samples.select_entries(150000,'random')
 assert len(a)==len(np.unique(a))==100000
 assert np.array_equal(a,samples.select_entries(150000,'random'))
 assert a.min()<50000 and a.max()>100000
 with pytest.raises(ValueError):samples.select_entries(10,'score_selected')

def test_cluster_uncertainty_keeps_pair_dependence():
 m=np.geomspace(2,120,20);mu=np.full(20,100.);mask=np.ones(20,bool)
 b=np.array([3,3,3,3])
 dependent=sweep.cluster_se(m,mu,mask,b,np.zeros(4))
 independent=sweep.cluster_se(m,mu,mask,b,np.arange(4))
 assert dependent==pytest.approx(2*independent)

def test_low_count_rows_are_preserved_for_every_model():
 d={'mass64':np.array([10.,15.]),'entry':np.array([1,2])}
 rows,_=sweep.evaluate(d,np.array([True,True]),'mass64','nominal','fixture','early','low')
 assert len(rows)==3
 assert all(r['status']=='LOW_COUNT' and r['in_range_pairs']==2 for r in rows)

def test_historical_and_sweep_detrending_match_on_synthetic_spectrum():
 from cms_wct.background_kill import evaluate_candidate,PairKillConfig
 from cms_wct.background_cv import BackgroundCandidate
 edges=np.geomspace(2,120,351);m=np.sqrt(edges[1:]*edges[:-1]);y=200*np.exp(-np.log(m)/2)*(1+.05*np.sin(np.log(m)))
 good=sweep.eligible(m)
 for degree in (7,12):
  _,expected=evaluate_candidate(m,y,BackgroundCandidate(f'cheb_d{degree}','chebyshev',degree),PairKillConfig(sweep.OMEGA,sweep.PHASE,sweep.MASKS))
  actual=sweep.fit(f'cheb_d{degree}',m,y,np.diff(edges),good)
  np.testing.assert_allclose(actual,expected,rtol=1e-10)
