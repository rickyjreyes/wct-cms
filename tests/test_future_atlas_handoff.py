"""Do not promote unperformed ATLAS development into a final signal."""
import json
import sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from screen_future_atlas_handoff import validate, ROOT, STATUS, main


def record():
    return json.loads(STATUS.read_text())


def candidate():
    return dict(repository="rickyjreyes/wct-atlas",
        source_commit="b"*40,freeze_commit="a"*40,
        manifest_sha256="c"*64,result_sha256="d"*64,
        coordinate="ln(m_mumu/1 GeV)",omega=7.025825825825827,
        prospectively_frozen=True,selection_mapping="audited independent map",
        detector_independence_assessment="separately reviewed",
        cms_frequency_statistic_evaluated=True,final_holdout_accessed=True,
        scientific_verdict="NEGATIVE",statistic=0.0,p_value=0.7)


def test_current_route_closed_is_not_signal():
    r=validate(record(),None)
    assert r["status"]=="EXTERNALLY_BLOCKED_FINAL_ATLAS_TEST_NOT_RUN"
    assert not r["new_atlas_signal_measured"]
    assert r["combined_significance"] is None


@pytest.mark.parametrize("mutation", [
    dict(cms_frequency_statistic_evaluated=False),
    dict(prospectively_frozen=False),
    dict(final_holdout_accessed=False),
    dict(coordinate="ln(q2)"),
    dict(omega=9.7),
    dict(manifest_sha256="fake"),
    dict(scientific_verdict="PASS"),
    dict(selection_mapping=""),
])
def test_unverified_or_wrong_claim_fails(mutation):
    v=candidate();v.update(mutation)
    with pytest.raises(ValueError):
        validate(record(),v)


def test_even_complete_candidate_is_only_unverified_screen(tmp_path):
    r=validate(record(),candidate())
    assert r["status"]=="CANDIDATE_REQUIRES_EXTERNAL_PROVENANCE_AND_PROTOCOL_REVIEW"
    assert r["claimed_verdict_unverified"]=="NEGATIVE"
    assert r["combined_significance"] is None
    assert r["external_hash_and_chronology_verified_by_this_script"] is False
    out=tmp_path/"screen.json"
    assert main(["--status",str(STATUS),"--out",str(out)])==0
    assert json.loads(out.read_text())["status"]=="EXTERNALLY_BLOCKED_FINAL_ATLAS_TEST_NOT_RUN"
    with pytest.raises(SystemExit):
        main(["--status",str(STATUS),"--out",str(out)])
