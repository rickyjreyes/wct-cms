#!/usr/bin/env python3
"""Local read-only ROOT-to-archive stage probe for CMS issue #9.

Runs only frozen extraction, histogram, continuum and locked-score diagnostics.
It does NOT re-run random nulls, replace archived evidence, or claim discovery.
Require a complete ROOT file verified against archived official CERN metadata.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import sys
import zlib
from pathlib import Path

# Executing `python scripts/trace_cms_verified_root_stages.py` makes scripts/
# sys.path[0]. The checkout uses a src-layout; pytest adds src automatically,
# but the standalone CLI must work even before `pip install -e .`.
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np
from cms_wct.cli import parser as cms_parser
from cms_wct.cmsio import extract_dimuon_masses, load_golden_json
from cms_wct.background import make_histogram, robust_smooth_background
from cms_wct.locked import fit_phase_locked_waveform
PREFIX = {"H1": ("30555", "127C"), "H2": ("30555", "183B"),
          "G1": ("30522", "05DD"), "G2": ("30522", "209D")}
MASKS = ("2.9:3.3", "3.55:3.85", "8.5:11.5", "80:100")
OMEGA = 7.025825825825827
PHASE = -0.2313916852932179


def published_identity(root: Path, label: str) -> dict:
    record, prefix = PREFIX[label]
    j = json.loads((root / f"audit/2026-09-06/baseline/cms{record}.json").read_text())
    matches = [f for f in j["metadata"]["_file_indices"][0]["files"]
               if f["filename"].startswith(prefix)]
    if len(matches) != 1:
        raise ValueError("Expected exactly one official file for frozen prefix")
    return matches[0]


def hash_input(path: Path, expected: dict) -> dict:
    if path.name != expected["filename"] or path.stat().st_size != expected["size"]:
        raise ValueError("Input filename/size does not match full official CERN ROOT; sparse caches rejected")
    sha = hashlib.sha256()
    checksum = 1
    with path.open("rb") as fh:
        while True:
            b = fh.read(8 * 1024 * 1024)
            if not b:
                break
            sha.update(b)
            checksum = zlib.adler32(b, checksum)
    adler = f"adler32:{checksum:08x}"
    if adler.lower() != expected["checksum"].lower():
        raise ValueError("Official full-file Adler-32 mismatch")
    return {"filename": path.name, "bytes": path.stat().st_size,
            "sha256": sha.hexdigest(), "adler32": adler,
            "official_checksum_verified": True}


def compare_arrays(reference, observed) -> dict:
    a, b = np.asarray(reference), np.asarray(observed)
    same_shape = a.shape == b.shape
    exact = bool(same_shape and np.array_equal(a,b, equal_nan=True))
    result = {"reference_shape": list(a.shape), "observed_shape": list(b.shape),
              "exact": exact, "max_abs_delta": None, "first_differing_index": None}
    if same_shape and a.size:
        diff = np.flatnonzero(~np.isclose(a.ravel(),b.ravel(),rtol=0,atol=0,equal_nan=True))
        result["first_differing_index"] = int(diff[0]) if len(diff) else None
        if np.issubdtype(a.dtype, np.number) and np.issubdtype(b.dtype, np.number):
            result["max_abs_delta"] = float(np.max(np.abs(a.astype(float)-b.astype(float))))
    return result


def stage_from_masses(masses, *, reference, root: Path):
    counts, edges, centers = make_histogram(masses,2,120,350,True)
    windows = [(float(s.split(":")[0]),float(s.split(":")[1])) for s in MASKS]
    model,bg = robust_smooth_background(centers,counts,7,6,3.5,windows)
    residual = (counts-model)/np.sqrt(np.maximum(model,1.0))
    mask=np.isfinite(residual) & (model>=5) & (centers>0)
    for lo,hi in windows: mask &= ~((centers>=lo)&(centers<=hi))
    x=np.log(centers[mask]); y=residual[mask]
    locked=fit_phase_locked_waveform(x,y,OMEGA,PHASE,positive_only=True)
    checks={}
    for name,a,b in (
        ("mass_bin_centers",reference["mass_center"],centers),
        ("histogram_counts",reference["count"],counts),
        ("background_fit_mask",reference["background_fit_mask"].astype(bool),bg),
        ("analysis_mask",reference["analysis_mask"].astype(bool),mask),
        ("background_mean",reference["background"],model),
        ("pearson_residuals",reference["residual"],residual)):
        checks[name]=compare_arrays(a,b)
    return {
        "checks_against_archived_audit_spectrum": checks,
        "observed_histogram_counts": [int(n) for n in counts],
        "analysis_bins":int(mask.sum()),
        "residual_rms":float(np.sqrt(np.mean(y*y))),
        "locked_signed_amplitude":float(locked.signed_amplitude),
        "locked_q":float(locked.delta_chi2)
    }


def diagnostic(root_file: Path, label: str, root: Path, step_sizes: list[str]) -> dict:
    official=published_identity(root,label)
    receipt=hash_input(root_file,official)
    cert=root/"data/cert/14221/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON_MuonPhys.txt"
    if not cert.is_file():
        raise ValueError("Frozen certification JSON missing")
    archive=root/f"audit/2026-09-06/baseline/{label}_masses.npz"
    spectrum=root/f"audit/2026-09-06/baseline/{label}/spectrum.csv"
    with np.load(archive,allow_pickle=False) as obj:
        audit_masses=obj["masses"].copy()
    audit_spectrum=np.genfromtxt(spectrum,delimiter=",",names=True)
    observed=[]
    for step in step_sizes:
        args=cms_parser().parse_args(["--input",str(root_file),
            "--golden-json",str(cert),"--mass-min","2","--mass-max","120",
            "--bins","350","--log-bins","--tight-id","--max-events","100000",
            "--step-size",step,"--fit-degree","7"])
        masses,counters=extract_dimuon_masses([str(root_file)],args,load_golden_json(str(cert)))
        observed.append({"step_size":step,"event_counts":counters,
           "mass_vector_vs_archived_audit":compare_arrays(audit_masses,masses),
           "pipeline":stage_from_masses(masses,reference=audit_spectrum,root=root),
           "_masses":masses})
    if len(observed)>1:
        pairwise=compare_arrays(observed[0]["_masses"],observed[1]["_masses"])
    else:
        pairwise=None
    for item in observed: del item["_masses"]
    same_vectors=all(v["mass_vector_vs_archived_audit"]["exact"] for v in observed)
    return {
      "schema":"CMS_VERIFIED_FULL_ROOT_STAGE_PROBE_V1",
      "label":label,"official_input_receipt":receipt,
      "certification_json_sha256":hashlib.sha256(cert.read_bytes()).hexdigest(),
      "archived_mass_npz_sha256":hashlib.sha256(archive.read_bytes()).hexdigest(),
      "archived_spectrum_sha256":hashlib.sha256(spectrum.read_bytes()).hexdigest(),
      "runs":observed,"between_step_sizes_mass_vector":pairwise,
      "all_step_sizes_match_archived_audit_mass_vector_exactly":same_vectors,
      "historical_full_intermediate_vectors_present":False,
      "cause_proven":False,
      "interpretation":("This isolates ROOT extraction and saved audit-stage differences. "
       "Agreement with an audit cache does not reconstruct unavailable historical "
       "full intermediate arrays; the cause must be supported by the observed "
       "first divergence and matching original historical provenance. No nulls scored.")
    }


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root",type=Path,default=ROOT)
    p.add_argument("--label",choices=sorted(PREFIX),required=True)
    p.add_argument("--root-file",type=Path,required=True,help="Complete original official CERN ROOT file")
    p.add_argument("--step-size",action="append",dest="steps",default=None)
    p.add_argument("--output",type=Path,required=True,help="Fresh output; never overwrites")
    a=p.parse_args(argv)
    if a.output.exists() or a.output.resolve()==a.root_file.resolve():
        p.error("Refusing to replace existing output or source input")
    receipt=diagnostic(a.root_file,a.label,a.root,a.steps or ["1 MB","100 MB"])
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print("Verified input; stage probe written to",a.output)
    print("Audit mass vector exact at all chunk sizes:",
          receipt["all_step_sizes_match_archived_audit_mass_vector_exactly"])
    return 0


if __name__=="__main__":
    raise SystemExit(main())
