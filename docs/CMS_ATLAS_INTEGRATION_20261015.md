# CMS–ATLAS dependency integration (#10; 2026-10-15)

**Classification: NO_VALID_ATLAS_SIGNAL_RESULT_TO_INTEGRATE.** This is
cross-repository provenance, NOT a rerun of ATLAS data and not a new CMS
replication. The machine-readable companion is
[CMS_ATLAS_INTEGRATION_STATUS_20261015.json](CMS_ATLAS_INTEGRATION_STATUS_20261015.json).

## Primary immutable records and preservation limits

- ATLAS exact V6.1 execution commit:
  [d5b3d7a3849453450a430c5b0ea46eaf456ed29f](https://github.com/rickyjreyes/wct-atlas/commit/d5b3d7a3849453450a430c5b0ea46eaf456ed29f)
- Preserved closure report and decision:
  [6511edc39308d588031e13081310575645573a3c](https://github.com/rickyjreyes/wct-atlas/commit/6511edc39308d588031e13081310575645573a3c),
  `audit/v61_closure_2026-10-04/decision.json`, ATLAS [#6](https://github.com/rickyjreyes/wct-atlas/issues/6),
  [#10](https://github.com/rickyjreyes/wct-atlas/issues/10).
- Input bundle SHA256:
  `c330891d79689fed6e85011b9725258c5ba09fe0fb3abc40bb813a8f6a6d645b`;
  bundle receipt SHA256:
  `fde752cef524d080a591e4d331aff02a7db6e6b1e22a7c4322043e2b260f9f4a`.
  Preserved execution artifact-manifest SHA256:
  `538d9df756764932e8b37afbd027a4485f2fea59bb842eb0c95970a8a3d08d0d`.
  These are recorded in the closure packet, not independently rehashed by
  this integration. Original artifact transfer remained incomplete (#6).
- Frozen development accounting: 256 allocated / 254 active / 2
  inactive zero-exposure. Candidate k=1..6; qualifying set `[]`.
  Classification **V6_FAILED_HIERARCHICAL_ROUTE_CLOSED**.
  **CMS-frequency target statistic NOT EVALUATED. Final ATLAS holdout NOT
  OPENED.** This is a negative *background-route eligibility* result, not
  evidence of zero CMS amplitude in ATLAS.
- Active-domain follow-on: frozen low-compute audit
  [af6da1f5cfdc0a6a685c79309db342dd2338aae9](https://github.com/rickyjreyes/wct-atlas/commit/af6da1f5cfdc0a6a685c79309db342dd2338aae9)
  and [#13–#15](https://github.com/rickyjreyes/wct-atlas/issues/15):
  `ACTIVE_DOMAIN_PROTOCOL_NOT_IDENTIFIABLE`,
  `NO_GO_FOR_NEW_ATLAS_RUN`, zero eligible cross-domain HE pairs.
  This did not evaluate a target or establish H0, HR, or HE.

## Observable, coordinate, and selection scope

The frozen CMS candidate `omega_CMS = 7.025825825825827` is in
`ell=ln(m_mumu/1 GeV)`, using opposite-sign muon pairs, pT≥4 GeV,
|eta|≤2.4, tight CMS ID, 2–120 GeV and 350 log bins with the declared
resonance/Z exclusions. The proposed ATLAS PHYSLITE mass-coordinate transfer
uses the **same** numerical omega only if ATLAS uses the same log-mass
coordinate. In a squared-mass coordinate `ln(q²)`, it becomes
`k=omega_CMS/2=3.512912912912913`. This is an identity of coordinates,
not empirical confirmation.

The ATLAS development builder did not enforce a one-for-one CMS Tight-muon
quality filter; trigger acceptance, reconstruction, object definitions,
weights and calibration also differ. CMS and ATLAS are different detector
chains, but **selection equivalence has not been established**, and the
absence of any valid final signal statistic precludes an independent-detector
replication claim. The developmental hierarchical background used frozen
stability/residual criteria and failed their joint adequacy. A later
active-domain hypothesis must not retroactively change that historical
rigid-transfer test.

| Hypothesis | What is on record | Permitted interpretation |
| --- | --- | --- |
| H0: no relation beyond background/null | Different domain-specific controls; no common calibrated comparison | Not declared proven |
| HR: rigid `k_W=k_0` | Frozen rigid ATLAS route closed in development; final signal untested | Route closed, **not** a measured negative ATLAS signal |
| HE: `kbar_W L_W=2 pi n_W` | Missing independent physical endpoint/winding/measurement rule; 0 admissible pairs | Protocol not identifiable; no ATLAS run authorized |

The CMS historical G2 fixed-waveform result remains reproducible *under its
declared pipeline* and Category C / significance background-dependent. Neither
an ATLAS route failure nor a later HE reformulation upgrades or downgrades
that CMS measurement by itself. No cross-detector p-value combination is
legitimate. See CMS #2, #4–#6 for the independent-background, efficiency and
simulation prerequisites. Keep CMS #10 tracked as externally blocked until a
new scientifically valid, preregistered ATLAS test with provenance is
actually completed.


## Successor final-result intake gate

The current V6.1 integration above is complete as documentation; the
original issue's **future final scientific result** is still externally
blocked. A separate non-scoring `scripts/screen_future_atlas_handoff.py`
provides a strict intake screen. It does not reinterpret the historical
V6.1 decision, execute ATLAS, verify source Git hashes over the network,
or authorize a new run.

```bash
python scripts/screen_future_atlas_handoff.py \
 --out results/atlas_handoff_not_run.json
# EXTERNALLY_BLOCKED_FINAL_ATLAS_TEST_NOT_RUN
```

Once a genuinely new, independently preregistered and authorized ATLAS final
result exists, the screen accepts a **copy** of a proposed machine-readable
handoff JSON via `--candidate ... --out NEW_OUTPUT`. Mandatory fields:
`repository`, full `source_commit`, full `freeze_commit`,
`manifest_sha256`, `result_sha256`, `coordinate`, `omega`,
`prospectively_frozen`, `selection_mapping`,
`detector_independence_assessment`, `cms_frequency_statistic_evaluated`,
`final_holdout_accessed`, `scientific_verdict` (POSITIVE, NEGATIVE or
INCONCLUSIVE), `statistic` and `p_value` (nullable only if not reported).
It refuses incorrect coordinate/frequency, missing source provenance or
attempted promotion of an unperformed test. Passing the schema yields ONLY
`CANDIDATE_REQUIRES_EXTERNAL_PROVENANCE_AND_PROTOCOL_REVIEW`, never
automatic acceptance, significance combination or README alteration.
The final result must be audited independently, and negative/inconclusive
outcomes retained. Original final-holdout content remains untouched while
the ATLAS protocol has a NO-GO gate.
