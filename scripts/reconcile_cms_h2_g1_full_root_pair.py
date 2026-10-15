#!/usr/bin/env python3
"""Recompute historical H2/G1 common-background CV from two full-ROOT probe receipts.

Requires new probe reports with exported per-step observed histogram counts.
Does not reconstruct unavailable historical intermediate vectors, calculate
new p-values, rewrite archives, or change the scientific Category C verdict.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
import numpy as np
from cms_wct.background import make_histogram
from cms_wct.background_kill import evaluate_pair_pipeline,PairKillConfig

WINDOWS=((2.9,3.3),(3.55,3.85),(8.5,11.5),(80.,100.))
CONFIG=PairKillConfig(omega=7.025825825825827,phase=-0.2313916852932179,
                      excluded_windows=WINDOWS)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validated_counts(payload, *, label, repo):
    if payload.get("schema")!="CMS_VERIFIED_FULL_ROOT_STAGE_PROBE_V1" or payload.get("label")!=label:
        raise ValueError("Wrong full-ROOT diagnostic schema/label")
    if payload.get("official_input_receipt",{}).get("official_checksum_verified") is not True:
        raise ValueError("Full official CERN input receipt missing")
    checks=payload.get("between_step_sizes_mass_vector")
    if not checks or checks.get("exact") is not True:
        raise ValueError("Ordered mass vectors must agree across full-source step sizes")
    runs=payload.get("runs",[])
    if {x.get("step_size") for x in runs}!={"1 MB","100 MB"} or len(runs)!=2:
        raise ValueError("Require exactly 1 MB and 100 MB full-source runs")
    cert=repo/"data/cert/14221/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON_MuonPhys.txt"
    if payload.get("certification_json_sha256")!=sha(cert):
        raise ValueError("Certification JSON mismatch")
    npz=repo/f"audit/2026-09-06/baseline/{label}_masses.npz"
    spectrum=repo/f"audit/2026-09-06/baseline/{label}/spectrum.csv"
    if payload.get("archived_mass_npz_sha256")!=sha(npz) or payload.get("archived_spectrum_sha256")!=sha(spectrum):
        raise ValueError("Archived audit source hashes mismatch")
    counts=[]
    for r in runs:
        c=np.asarray(r.get("pipeline",{}).get("observed_histogram_counts"),dtype=float)
        if c.shape!=(350,) or np.any(~np.isfinite(c)) or np.any(c<0) or np.any(c!=np.floor(c)):
            raise ValueError("Invalid or absent count-vector export; rerun upgraded full-ROOT probe")
        counts.append(c)
    if not np.array_equal(counts[0],counts[1]):
        raise ValueError("Histogram counts disagree across step sizes")
    return counts[0]


def replay(h2, g1, *, repo=ROOT):
    a=validated_counts(h2,label="H2",repo=repo)
    b=validated_counts(g1,label="G1",repo=repo)
    # Use the exact nominal histogram-center implementation, independent of archive serialization.
    _,_,centers=make_histogram(np.empty(0),2,120,350,True)
    observed,_,_=evaluate_pair_pipeline(centers,a,centers,b,None,CONFIG)
    archived=json.loads((repo/"audit/2026-09-06/baseline/pair_observed.json").read_text())
    archived_a=np.genfromtxt(repo/"audit/2026-09-06/baseline/H2/spectrum.csv",delimiter=",",names=True)
    archived_b=np.genfromtxt(repo/"audit/2026-09-06/baseline/G1/spectrum.csv",delimiter=",",names=True)
    audit_replay,_,_=evaluate_pair_pipeline(centers,archived_a["count"],centers,archived_b["count"],None,CONFIG)
    return dict(
      schema="CMS_H2_G1_FULL_SOURCE_PAIR_RECONCILIATION_V1",
      frozen_config=CONFIG.to_dict(),
      verified_h2_root_receipt=h2["official_input_receipt"],
      verified_g1_root_receipt=g1["official_input_receipt"],
      full_source_pair=observed,
      archived_audit_pair_score=float(archived["primary_pair_score_min_locked_delta_chi2"]),
      archived_audit_selected_background=archived["selected_background"]["name"],
      audit_counts_recomputed_pair_score=float(audit_replay["primary_pair_score_min_locked_delta_chi2"]),
      full_source_vs_archived_h2_histogram_equal=bool(np.array_equal(a,archived_a["count"])),
      full_source_vs_archived_g1_histogram_equal=bool(np.array_equal(b,archived_b["count"])),
      historic_pair_q_approx_from_existing_docs=108.49776,
      historic_pair_q_is_approximate=True,
      historical_pair_full_intermediate_arrays_available=False,
      cause_proven=False,
      scientific_status="CATEGORY_C_BACKGROUND_DEPENDENT_UNCHANGED",
      no_new_nulls_or_p_values=True,
    )


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--h2-probe",type=Path,required=True)
    p.add_argument("--g1-probe",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--root",type=Path,default=ROOT)
    a=p.parse_args(argv)
    if a.output.exists() or a.output.resolve() in (a.h2_probe.resolve(),a.g1_probe.resolve()):
        p.error("Refuse overwrite of any existing evidence/input")
    h=json.loads(a.h2_probe.read_text())
    g=json.loads(a.g1_probe.read_text())
    report=replay(h,g,repo=a.root)
    report["h2_probe_sha256"]=sha(a.h2_probe)
    report["g1_probe_sha256"]=sha(a.g1_probe)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("Paired-CV result (fresh full-source counts):",
          report["full_source_pair"]["primary_pair_score_min_locked_delta_chi2"])
    print("Selected background:",report["full_source_pair"]["selected_background"]["name"])
    print("Archived audit score:",report["archived_audit_pair_score"])
    print("Cause proven:",report["cause_proven"])
    return 0


if __name__=="__main__":
    raise SystemExit(main())
