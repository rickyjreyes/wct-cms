# Execute the remaining CMS issue #9 numerical cause test (2026-10-15)

The archived summary-only comparator is already implemented and passed locally
on the BLD machine (38/38 after using a writable pytest base-temp). It finds
the first **saved** G2 discrepancy in residual RMS. This **does not** prove
the first arithmetic/source difference.

This new probe runs the next genuinely discriminating test on **the full
official CERN ROOT input** and checks both archived official size and
Adler-32 before reading events. It rejects partial/sparse-cache files. The
complete G2 official file is:

`209D94D9-B6D5-A34B-A2A3-CBB7E4EA8ADF.root`
- CERN record `30522`, 2,284,233,920 bytes, Adler-32 `b18e8b13`.
- Frozen golden JSON, first 100,000 input entries before certification,
  Tight-muon pair cuts, 2–120 GeV, 350 log bins, frozen resonance exclusions,
  degree-seven Chebyshev and omega/phase are reused without retuning.
- The original archived audit `G2_masses.npz` and `G2/spectrum.csv` are
  *read only*. Full ROOT SHA256 and golden JSON SHA256 are emitted in a
  **new** report. No null p-values or physics discoveries are produced.

Use Git Bash, adjusting the source path to your **complete existing** official
file (never supply `.partial.root`). Results/ is untracked output; use a
fresh filename rather than overwriting historical artifacts:

```bash
cd /r/RESEARCH/wct-cms
python scripts/trace_cms_verified_root_stages.py \
  --label G2 \
  --root-file 'R:/RESEARCH/wct-cms/data/audit_inputs/209D94D9-B6D5-A34B-A2A3-CBB7E4EA8ADF.root' \
  --output results/cms_g2_verified_root_stage_probe.json
python -m pytest -q --basetemp='R:/RESEARCH/.pytest_wct_cms_tmp'
```

If no complete file is available, obtain it through the existing
`experiments/cms_audit/fetch_inputs.py` only with sufficient disk space;
that downloader fetches four large files and verifies official Adler-32.
Or use your already-downloaded complete official G2 file wherever it resides.

Repeat for H2 and G1 if necessary with labels and matching filenames from
`published_identity`. The probe compares selected **float64 mass vectors in
order**, histogram bin assignments, background-fit and analysis masks,
model predictions, Pearson residuals, and frozen score for both `1 MB` and
`100 MB` uproot step sizes. A first mismatch in mass-vector order, counts,
mask, or continuum provides a *candidate first reproducible computational
difference*; evaluate its actual cause against historical dependencies and
original intermediate arrays before closing #9. Equal event totals alone
never certify identical mass arrays. If audit and both full-file steps agree
exactly but the historical score still differs, prioritize dependency and
historical-fit convergence analysis rather than blaming step size.

Historical source/intermediate evidence and the modern Category C model
dependence are independent of this numerical debugging. Preserve negative and
inconclusive comparisons; never overwrite historical arrays to force equality.

## 2026-10-15 BLD outcome and next executed gates

User independently ran G2 against complete official full ROOT:
`209D94D9-B6D5-A34B-A2A3-CBB7E4EA8ADF.root`,
2,284,233,920 bytes, official Adler-32 `b18e8b13`,
full-file SHA256 `91b22c5897c0b4ff3d50bf6dd1fb21a3d527fa8ac5bdfbd74e6f01bd7ad401f1`.
Both 1MB and 100MB selected mass vectors were exactly equal in default BLD
and the isolated pinned Python 3.12 audit environment. All four yielded
historical G2 residual RMS `1.729361458528311`, locked amplitude
`0.970861774600961`, Q `126.2832399542242`, not the archived
audit Q `126.27500073555109`. The archived G2 mass NPZ **exactly**
reconstructs the archived spectrum's 350 histogram counts and 38,917 in-range
masses. Fresh vs archived NPZ has 4,547/45,803 differing ordered masses,
maximum 0.025537729263305664 GeV; sorted vectors also differ. The first saved
divergence is therefore the mass vector, not isolated histogram or fitting
logic. Neither the old sparse-cache bytes nor old event keys have been
revalidated. Do not assert sparse-cache corruption or a unique cause.

### A. Controlled G2 arithmetic precision diagnostic (new)

Use the pinned audit environment already created locally, or another environment
with `pip install -e .`; use a fresh output path:

```bash
OPENBLAS_NUM_THREADS=1 python scripts/probe_cms_g2_mass_precision.py \
  --root-file 'R:/RESEARCH/wct-cms/data/audit_inputs/209D94D9-B6D5-A34B-A2A3-CBB7E4EA8ADF.root' \
  --output results/cms_g2_same_pair_precision_20261015.json
```

The script first verifies the *entire* file against official CERN size and
Adler-32. It reads exactly the first 100,000 original Events entries and
applies the frozen golden JSON, Tight ID, pT, eta and all OS-pair selection
once, comparing nominal branch-dtype arithmetic to explicit float64 promotion
of all four muon-vector **components before arithmetic** on the same pairs.
It reports ordered/sorted discrepancies, first differing masses, quantiles and
float32-exact fractions against the immutable NPZ. It is NOT a validation of
old sparse cache bytes; even a better match is a causal lead, not proof.

### B. H2/G1 full-source paired-score reconciliation (new)

The full-ROOT stage probe now includes `observed_histogram_counts` in each
run to preserve actual source-verified 350-bin counts. Regenerate **fresh**
H2 and G1 probes with their complete official original files, e.g.:

```bash
OPENBLAS_NUM_THREADS=1 python scripts/trace_cms_verified_root_stages.py \
  --label H2 --root-file 'R:/RESEARCH/wct-cms/data/audit_inputs/183BFB78-7B5E-734F-BBF5-174A73020F89.root' \
  --output results/cms_h2_full_root_pair_probe.json
OPENBLAS_NUM_THREADS=1 python scripts/trace_cms_verified_root_stages.py \
  --label G1 --root-file 'R:/RESEARCH/wct-cms/data/audit_inputs/05DD095C-F6C3-9A4F-9FB3-348A5A6403D5.root' \
  --output results/cms_g1_full_root_pair_probe.json

OPENBLAS_NUM_THREADS=1 python scripts/reconcile_cms_h2_g1_full_root_pair.py \
  --h2-probe results/cms_h2_full_root_pair_probe.json \
  --g1-probe results/cms_g1_full_root_pair_probe.json \
  --output results/cms_h2_g1_full_source_pair_reconciliation.json
```

The replay requires matching certification/audit file SHA256 receipts, official
checksum-verified provenance, 1MB/100MB exact ordered-mass and count parity,
and independently recomputes the **original WCT-blind common-background
15-candidate blocked-CV selection** and frozen paired statistic from new
observed counts. It compares with the archived audit paired result and also
reruns the archived **counts** through identical paired CV; the approximate
historical pair Q=108.49776 is marked approximate, not a fabricated exact
archive intermediate. It performs no pseudoexperiment, no p-value combining,
no output overwrite and no physics reinterpretation. If either full official
file is unavailable, record **EXTERNALLY_BLOCKED_INPUT_UNAVAILABLE** locally
and do not use sparse file caches as substitutes for checksum verification.

### Closure boundary

#9 remains OPEN until same-pair precision and/or original-cache provenance
sufficiently bounds the G2 stage cause, the actual H2/G1 paired-score
discrepancy is reconstructed under verified full-file inputs, and an
appropriate numerical regression and scientific-status annotation are
archived. The larger conventional-background identifiability issue is
independent: CMS stays Category C regardless of a numerical fix.
