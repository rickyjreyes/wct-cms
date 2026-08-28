from types import SimpleNamespace

import awkward as ak

from cms_wct import cmsio


def test_tight_id_mask_composes_without_inplace_awkward_ufunc(monkeypatch):
    arrays = {
        "run": ak.Array([1, 1]),
        "luminosityBlock": ak.Array([1, 1]),
        "Muon_pt": ak.Array([[10.0, 12.0], [8.0, 9.0]]),
        "Muon_eta": ak.Array([[0.1, -0.2], [0.3, -0.4]]),
        "Muon_phi": ak.Array([[0.0, 3.141592653589793], [0.2, 2.8]]),
        "Muon_mass": ak.Array([[0.105, 0.105], [0.105, 0.105]]),
        "Muon_charge": ak.Array([[1, -1], [1, -1]]),
        "Muon_tightId": ak.Array([[True, True], [True, False]]),
    }

    def fake_iterate(*args, **kwargs):
        yield arrays

    monkeypatch.setattr(cmsio.uproot, "iterate", fake_iterate)

    args = SimpleNamespace(
        tree="Events",
        step_size="100 MB",
        max_events=None,
        muon_pt_min=4.0,
        muon_eta_max=2.4,
        tight_id=True,
        medium_id=False,
    )

    masses, counters = cmsio.extract_dimuon_masses(["dummy.root"], args)

    assert len(masses) == 1
    assert counters["events_read"] == 2
    assert counters["events_after_json"] == 2
    assert counters["selected_muons"] == 3
    assert counters["opposite_sign_pairs"] == 1
