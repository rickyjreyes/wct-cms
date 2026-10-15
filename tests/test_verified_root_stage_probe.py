"""No-root CI coverage for the full-file local causal probe."""
import hashlib
import sys
import zlib
from pathlib import Path
import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from trace_cms_verified_root_stages import compare_arrays, hash_input, published_identity, ROOT


def test_official_identity_for_g2():
    identity=published_identity(ROOT,"G2")
    assert identity["filename"]=="209D94D9-B6D5-A34B-A2A3-CBB7E4EA8ADF.root"
    assert identity["size"]==2284233920
    assert identity["checksum"]=="adler32:b18e8b13"


def test_stage_comparison_finds_first_mass_difference():
    a=np.array([2.,3.,4.]); b=np.array([2.,3.000001,4.])
    out=compare_arrays(a,b)
    assert out["exact"] is False and out["first_differing_index"]==1
    assert out["max_abs_delta"]==pytest.approx(1e-6)
    assert compare_arrays(a,a)["exact"] is True
    assert compare_arrays(a,a[:2])["exact"] is False


def test_official_size_name_and_adler_guard(tmp_path):
    raw=b"full synthetic fixture, not a CMS event file"
    file=tmp_path/"synthetic.root"; file.write_bytes(raw)
    metadata=dict(filename=file.name,size=len(raw),
                  checksum=f"adler32:{zlib.adler32(raw):08x}")
    r=hash_input(file,metadata)
    assert r["official_checksum_verified"] is True
    assert r["sha256"]==hashlib.sha256(raw).hexdigest()
    with pytest.raises(ValueError,match="filename/size"):
        hash_input(file,{**metadata,"size":len(raw)+1})
    with pytest.raises(ValueError,match="Adler"):
        hash_input(file,{**metadata,"checksum":"adler32:00000000"})


def test_legacy_counts_dont_prove_mass_vector_equality():
    assert compare_arrays([2.,3.],[2.1,2.9])["exact"] is False


def test_probe_standalone_help_from_clean_checkout():
    """Direct invocation must work without pytest's implicit src path injection."""
    import os
    import subprocess
    env=os.environ.copy()
    env.pop("PYTHONPATH",None)
    # -I ignores user site / PYTHONPATH; the probe must inject local src itself.
    run=subprocess.run([sys.executable,"-I",
        str(ROOT/"scripts/trace_cms_verified_root_stages.py"),"--help"],
        cwd=str(ROOT),env=env,capture_output=True,text=True,timeout=30)
    assert run.returncode==0,run.stderr
    assert "--root-file" in run.stdout
