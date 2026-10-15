"""Reproduction guard: earliest differing archived field is not a proven cause."""
import json
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"scripts"))
from compare_cms_numeric_reproduction import ROOT, HIST_G2, AUDIT_G2, AUDIT_PAIR, compare, main


def test_stagewise_comparison_preserves_status():
    h=json.loads((ROOT/HIST_G2).read_text())
    a=json.loads((ROOT/AUDIT_G2).read_text())
    p=json.loads((ROOT/AUDIT_PAIR).read_text())
    result=compare(h,a,p)
    assert result["matched_selected_counts"] is True
    assert result["first_differing_archived_stage"]=="residual_summary"
    assert result["first_differing_archived_field"]=="residual_rms"
    assert result["observed_pair_winner"]==result["historical_pair_winner"]=="spline_s2"
    assert result["pair_relative_difference_pct"] < 0.06
    assert result["step_size_is_a_cause_proven"] is False
    assert result["first_computational_divergence_proven"] is False
    assert result["physical_significance_inferred"] is False


def test_readonly_cli_preserves_inputs(tmp_path):
    sources=(HIST_G2,AUDIT_G2,AUDIT_PAIR)
    before={s:(ROOT/s).read_bytes() for s in sources}
    out=tmp_path/"comparison.json"
    assert main(["--root",str(ROOT),"--out",str(out)])==0
    assert out.exists()
    assert {s:(ROOT/s).read_bytes() for s in sources}==before
