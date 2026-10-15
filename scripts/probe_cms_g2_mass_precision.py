#!/usr/bin/env python3
"""Read-only G2 controlled mass-precision replay for CMS numerical issue #9.

Use the verified complete CERN ROOT input and the same selected pair ordering.
Compare archived masses with nominal (branch-dtype arithmetic) and explicitly
float64-promoted *input components before arithmetic*. A better match is
diagnostic, not proof of old sparse-cache provenance or original machine mode.
No historical artifacts are altered and no nulls are evaluated.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import awkward as ak
import numpy as np
import uproot
from cms_wct.cmsio import golden_mask, invariant_mass_from_pairs, load_golden_json
from trace_cms_verified_root_stages import hash_input, published_identity, compare_arrays


def select_pairs(arrays, golden, *, promoted: bool):
    gm = ak.Array(golden_mask(ak.to_numpy(arrays["run"]),
                              ak.to_numpy(arrays["luminosityBlock"]), golden))
    pt = arrays["Muon_pt"][gm]; eta = arrays["Muon_eta"][gm]
    keep = (pt >= 4.) & (abs(eta) <= 2.4) & arrays["Muon_tightId"][gm]
    def value(branch):
        selected = arrays[branch][gm][keep]
        return ak.values_astype(selected, np.float64) if promoted else selected
    mu = ak.zip(dict(pt=value("Muon_pt"), eta=value("Muon_eta"),
                     phi=value("Muon_phi"), mass=value("Muon_mass"),
                     charge=arrays["Muon_charge"][gm][keep]))
    pair = ak.combinations(mu, 2, fields=["a", "b"])
    return pair[(pair.a.charge * pair.b.charge) < 0]


def mass_vector(pairs):
    v = ak.to_numpy(ak.flatten(invariant_mass_from_pairs(pairs), axis=None))
    return v[np.isfinite(v)].astype(np.float64, copy=False)


def exact_f32_fraction(v):
    a = np.asarray(v, dtype=np.float64)
    return int(np.count_nonzero(a == a.astype(np.float32).astype(np.float64)))


def summarize(archived, candidate):
    x = compare_arrays(archived, candidate)
    if archived.shape != candidate.shape:
        raise ValueError("Different vector shapes: selection/pair order must be reconciled first")
    change = archived != candidate
    delta = np.abs(archived - candidate)
    x.update(
        differing_ordered_elements=int(np.count_nonzero(change)),
        identical_ordered_elements=int(np.count_nonzero(~change)),
        same_sorted_values=bool(np.array_equal(np.sort(archived), np.sort(candidate))),
        max_sorted_abs_delta=float(np.max(np.abs(np.sort(archived)-np.sort(candidate)))),
        median_abs_delta=float(np.median(delta)),
        quantiles_abs_delta=[float(v) for v in np.quantile(delta,[.5,.9,.99,1.])],
        candidate_exact_float32_count=exact_f32_fraction(candidate),
        first_differences=[dict(index=int(i), archived=float(archived[i]),
                                candidate=float(candidate[i]),
                                delta=float(candidate[i]-archived[i]))
                           for i in np.flatnonzero(change)[:12]],
    )
    return x


def run_probe(root_file, repo):
    receipt = hash_input(root_file, published_identity(repo, "G2"))
    cert = repo / "data/cert/14221/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON_MuonPhys.txt"
    archive = repo / "audit/2026-09-06/baseline/G2_masses.npz"
    golden = load_golden_json(str(cert))
    branches = ["run", "luminosityBlock", "Muon_pt", "Muon_eta",
                "Muon_phi", "Muon_mass", "Muon_charge", "Muon_tightId"]
    with uproot.open(root_file) as source:
        arrays = source["Events"].arrays(branches, entry_start=0,
                                          entry_stop=100000, library="ak", how=dict)
    if len(arrays["run"]) != 100000:
        raise ValueError("The official full ROOT file does not expose 100000 first entries")
    nominal = mass_vector(select_pairs(arrays, golden, promoted=False))
    promoted = mass_vector(select_pairs(arrays, golden, promoted=True))
    with np.load(archive, allow_pickle=False) as z:
        archived = z["masses"].copy()
    if nominal.shape != archived.shape or promoted.shape != archived.shape:
        raise ValueError("Selected pair count/order requires further event-key investigation")
    return dict(
      schema="CMS_G2_SAME_PAIR_PRECISION_DIAGNOSTIC_V1",
      official_input_receipt=receipt,
      golden_json_sha256=hashlib.sha256(cert.read_bytes()).hexdigest(),
      archived_mass_npz_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),
      pair_count=int(len(nominal)),
      archived_exact_float32_count=exact_f32_fraction(archived),
      nominal_vs_archive=summarize(archived, nominal),
      promoted_float64_vs_archive=summarize(archived, promoted),
      promoted_vs_nominal=summarize(nominal,promoted),
      cause_proven=False,
      event_identity_independently_checked=False,
      physics_signal_scored=False,
      interpretation=("Branches and pair ordering are held fixed within the fresh run. "
        "This tests a component-promotion arithmetic hypothesis only. Even exact "
        "numeric agreement cannot prove how archived sparse-cache bytes were "
        "formed; event-level provenance and H2/G1 pair CV still require review."),
    )


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root-file",required=True,type=Path)
    p.add_argument("--root",type=Path,default=ROOT)
    p.add_argument("--output",required=True,type=Path)
    a=p.parse_args(argv)
    if a.output.exists() or a.output.resolve()==a.root_file.resolve():
        p.error("Refusing to overwrite any existing result or ROOT input")
    result=run_probe(a.root_file,a.root)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("Read-only precision report:",a.output)
    print("Nominal matches archived:",result["nominal_vs_archive"]["exact"])
    print("Promoted matches archived:",result["promoted_float64_vs_archive"]["exact"])
    print("Cause proven:",result["cause_proven"])
    return 0


if __name__=="__main__":
    raise SystemExit(main())
