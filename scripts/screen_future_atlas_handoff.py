#!/usr/bin/env python3
"""Issue #10: fail-closed non-scoring intake for a *future* final ATLAS result.

The original V6.1 route-closed decision is never reclassified. This script
does not run an ATLAS model or verify external commits/hashes on its own.
An acceptable input only becomes a candidate for independent evidence review.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STATUS = ROOT/"docs/CMS_ATLAS_INTEGRATION_STATUS_20261015.json"
SHA40 = re.compile(r"^[0-9a-f]{40}$")
SHA64 = re.compile(r"^[0-9a-f]{64}$")
VERDICTS = {"POSITIVE","NEGATIVE","INCONCLUSIVE"}


def validate(current: dict, candidate: dict | None) -> dict:
    atlas=current["atlas"]
    if atlas["cms_frequency_statistic_evaluated"] or atlas["final_holdout_accessed"] or (
        atlas["final_signal_result"] is not None):
        raise ValueError("Historical V6.1 record changed; refuse integration")
    if candidate is None:
        return {"schema":"CMS_ATLAS_FINAL_HANDOFF_SCREEN_V1",
                "status":"EXTERNALLY_BLOCKED_FINAL_ATLAS_TEST_NOT_RUN",
                "new_atlas_signal_measured":False,"historical_v61_retained":True,
                "combined_significance":None}
    required=("repository","source_commit","freeze_commit","manifest_sha256",
              "result_sha256","coordinate","omega","prospectively_frozen",
              "selection_mapping","detector_independence_assessment",
              "cms_frequency_statistic_evaluated","final_holdout_accessed",
              "scientific_verdict","statistic","p_value")
    missing=[k for k in required if k not in candidate]
    if missing:
        raise ValueError("Missing future result provenance: "+", ".join(missing))
    if candidate["repository"]!="rickyjreyes/wct-atlas":
        raise ValueError("Wrong independent source repository")
    for k in ("source_commit","freeze_commit"):
        if not isinstance(candidate[k],str) or not SHA40.fullmatch(candidate[k]):
            raise ValueError("Invalid immutable source/freeze commit")
    for k in ("manifest_sha256","result_sha256"):
        if not isinstance(candidate[k],str) or not SHA64.fullmatch(candidate[k]):
            raise ValueError("Invalid artifact digest")
    if candidate["source_commit"]==candidate["freeze_commit"]:
        raise ValueError("Cannot infer prospective chronology from one identical commit")
    if candidate["prospectively_frozen"] is not True or (
        candidate["cms_frequency_statistic_evaluated"] is not True or
        candidate["final_holdout_accessed"] is not True):
        raise ValueError("A development route or unperformed test is not a final signal result")
    if candidate["coordinate"]!="ln(m_mumu/1 GeV)" or not math.isclose(
          float(candidate["omega"]),float(current["cms"]["omega"]),rel_tol=0,abs_tol=1e-10):
        raise ValueError("Not the originally frozen rigid CMS mass-coordinate test")
    if candidate["scientific_verdict"] not in VERDICTS:
        raise ValueError("Must preserve positive, negative, or inconclusive outcome explicitly")
    if not isinstance(candidate["selection_mapping"],str) or not candidate["selection_mapping"].strip():
        raise ValueError("Missing explicit selection mapping")
    if not isinstance(candidate["detector_independence_assessment"],str) or not candidate["detector_independence_assessment"].strip():
        raise ValueError("Missing independence assessment")
    if candidate["statistic"] is not None and not math.isfinite(float(candidate["statistic"])):
        raise ValueError("Nonfinite statistic")
    if candidate["p_value"] is not None and not (0 <= float(candidate["p_value"]) <= 1):
        raise ValueError("Invalid reported p-value")
    return {"schema":"CMS_ATLAS_FINAL_HANDOFF_SCREEN_V1",
            "status":"CANDIDATE_REQUIRES_EXTERNAL_PROVENANCE_AND_PROTOCOL_REVIEW",
            "claimed_verdict_unverified":candidate["scientific_verdict"],
            "claimed_source_commit":candidate["source_commit"],
            "historical_v61_retained":True,
            "external_hash_and_chronology_verified_by_this_script":False,
            "independence_for_combination_established":False,
            "combined_significance":None,
            "next_action":"Independently check exact source/artifact bytes and timestamps, selection mapping, frozen protocol, adequate null/calibration, and whether a new test was actually authorized. Preserve negative and inconclusive outcomes; do not auto-update README or combine p-values."}


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--status",type=Path,default=STATUS)
    p.add_argument("--candidate",type=Path)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args(argv)
    if a.out.exists() or a.out.resolve() in {z.resolve() for z in (a.status,a.candidate) if z is not None}:
        p.error("Refusing to overwrite source or prior screen")
    source=json.loads(a.status.read_text())
    candidate=json.loads(a.candidate.read_text()) if a.candidate else None
    result=validate(source,candidate)
    result["frozen_cms_integration_sha256"]=hashlib.sha256(a.status.read_bytes()).hexdigest()
    if a.candidate:
        result["candidate_sha256"]=hashlib.sha256(a.candidate.read_bytes()).hexdigest()
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(result["status"],"->",a.out)
    return 0


if __name__=="__main__":
    raise SystemExit(main())
