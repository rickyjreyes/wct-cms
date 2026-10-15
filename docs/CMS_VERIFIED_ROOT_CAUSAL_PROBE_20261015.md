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
