"""Same-pair precision diagnostic: no external CERN file required for CI."""
import sys
from pathlib import Path
import numpy as np
import awkward as ak
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"scripts"))
from probe_cms_g2_mass_precision import select_pairs, mass_vector, summarize, exact_f32_fraction


def test_float32_exact_fraction_and_quantiles():
    archived=np.array([1.0, 2.0, 3.0],dtype=np.float64)
    fresh=np.array([1.0, 2.00001, 3.0],dtype=np.float64)
    a=summarize(archived,fresh)
    assert a["differing_ordered_elements"]==1
    assert a["identical_ordered_elements"]==2
    assert a["same_sorted_values"] is False
    assert a["first_differences"][0]["index"]==1
    assert exact_f32_fraction(np.array([1., 2., 0.1],float))==2


def test_same_selected_pairs_and_explicit_promotion():
    # Two muons in certified event, one in uncertified event.
    a={"run":ak.Array([1,2]),"luminosityBlock":ak.Array([1,1]),
       "Muon_pt":ak.values_astype(ak.Array([[4.5,4.1],[6.]]),np.float32),
       "Muon_eta":ak.values_astype(ak.Array([[.1,-.2],[.1]]),np.float32),
       "Muon_phi":ak.values_astype(ak.Array([[0.,3.1],[0.]]),np.float32),
       "Muon_mass":ak.values_astype(ak.Array([[.105,.105],[.105]]),np.float32),
       "Muon_charge":ak.Array([[1,-1],[1]]),
       "Muon_tightId":ak.Array([[True,True],[True]])}
    golden={1:[(1,1)]}
    nominal=select_pairs(a,golden,promoted=False)
    high=select_pairs(a,golden,promoted=True)
    assert len(mass_vector(nominal))==len(mass_vector(high))==1
    assert ak.to_numpy(ak.flatten(high.a.pt,axis=None)).dtype==np.float64
    assert ak.to_numpy(ak.flatten(nominal.a.pt,axis=None)).dtype==np.float32
    with pytest.raises(ValueError,match="Different vector shapes"):
        summarize(np.array([1.]),np.array([1.,2.]))
