"""No-ROOT replay guards for G2-independent H2/G1 pair result."""
import sys
from pathlib import Path
import json
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from reconcile_cms_h2_g1_full_root_pair import validated_counts, sha, ROOT


def stub(label):
    base=ROOT/"audit/2026-09-06/baseline"
    return {
        "schema":"CMS_VERIFIED_FULL_ROOT_STAGE_PROBE_V1",
        "label":label,
        "official_input_receipt":{"official_checksum_verified":True},
        "between_step_sizes_mass_vector":{"exact":True},
        "certification_json_sha256":sha(ROOT/"data/cert/14221/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON_MuonPhys.txt"),
        "archived_mass_npz_sha256":sha(base/f"{label}_masses.npz"),
        "archived_spectrum_sha256":sha(base/f"{label}/spectrum.csv"),
        "runs":[{"step_size":step,"pipeline":{"observed_histogram_counts":[0]*350}}
                for step in ("1 MB","100 MB")]
    }


def test_replay_rejects_unverified_or_changed_chunk_counts():
    s=stub("H2")
    assert validated_counts(s,label="H2",repo=ROOT).shape==(350,)
    s["runs"][1]["pipeline"]["observed_histogram_counts"][13]=1
    with pytest.raises(ValueError,match="disagree"):
        validated_counts(s,label="H2",repo=ROOT)
    s=stub("H2");s["official_input_receipt"]["official_checksum_verified"]=False
    with pytest.raises(ValueError,match="receipt"):
        validated_counts(s,label="H2",repo=ROOT)
    s=stub("H2");s["archived_mass_npz_sha256"]="altered"
    with pytest.raises(ValueError,match="source hashes"):
        validated_counts(s,label="H2",repo=ROOT)
    s=stub("H2");s["runs"][0]["pipeline"].pop("observed_histogram_counts")
    with pytest.raises(ValueError,match="absent"):
        validated_counts(s,label="H2",repo=ROOT)
