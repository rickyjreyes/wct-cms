"""Issue #10: route closure must never be serialized as an ATLAS signal test."""
import json
from pathlib import Path


def test_atlas_provenance_semantics():
    root=Path(__file__).resolve().parents[1]
    j=json.loads((root/"docs/CMS_ATLAS_INTEGRATION_STATUS_20261015.json").read_text())
    assert j["atlas"]["v61_classification"]=="V6_FAILED_HIERARCHICAL_ROUTE_CLOSED"
    assert j["atlas"]["qualifying_k"]==[]
    assert j["atlas"]["cms_frequency_statistic_evaluated"] is False
    assert j["atlas"]["final_holdout_accessed"] is False
    assert j["atlas"]["final_signal_result"] is None
    assert j["atlas"]["combined_p_value"] is None
    assert j["atlas"]["active_domain_protocol"]=="ACTIVE_DOMAIN_PROTOCOL_NOT_IDENTIFIABLE"
    assert j["atlas"]["active_domain_admissible_HE_pairs"]==0
    assert j["atlas"]["comparison"]["selection_equivalence_established"] is False
    assert j["atlas"]["comparison"]["squared_mass_coordinate_frequency"] == (
        j["cms"]["omega"]/2)
