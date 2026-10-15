# CMS numerical audit: archived-stage comparison (#9, 2026-10-15)

Run `python scripts/compare_cms_numeric_reproduction.py` from repository root.
Output is a **new** derived `results/cms_numeric_archive_comparison.json`, not a
replacement for original evidence. The report compares archived G2 historical
`docs/CMS_RUN2016G_FILE2_PHASE_LOCK_RESULT_2026-08-31.json` with independently
rerun audit `audit/2026-09-06/baseline/G2_locked/summary.json`, and the preserved
audit pair score. The subsequent external AI reviewer read these archived
outputs but did **not** rerun ROOT inputs.

| Comparable checkpoint | Historical | Audit | Reading |
| --- | ---: | ---: | --- |
| G2 events read/after certification | 100,000 / 100,000 | 100,000 / 100,000 | Match |
| selected muons / OS pairs | 143,052 / 45,803 | 143,052 / 45,803 | Match |
| masses in 2–120 GeV / analysis bins | 38,917 / 287 | 38,917 / 287 | Match |
| first differing **archived field**: residual RMS | 1.7293614585 | 1.7289454627 | Differ |
| locked amplitude | 0.9708617746 | 0.9708301027 | Differ slightly |
| locked Q | 126.2832399542 | 126.2750007356 | Differ slightly |
| selected H2/G1 pair Q | ≈108.49776 | 108.4374621327 | Difference ≈0.056% |
| pair-selected background | spline_s2 | spline_s2 | Match |

**This does not identify the first computational point of divergence.**
The saved results are summaries rather than a complete matched vector of
historical mass arrays, fit masks, intermediate baselines, and residuals.
Comparing counts alone cannot establish bitwise event-vector identity. The
archived historical config uses `100 MB` input steps, the sparse-cache audit
uses `1 MB`; this is an observed environment difference, **not a proven cause**.
The audit environment is Python 3.12.13, NumPy 2.3.5, SciPy 1.17.0, uproot
5.7.6, awkward 2.13.0; the original historical full dependency lock and
input checksums are not independently reconstructed here.

## Remaining causal experiment

1. Verify official file identities, original full ROOT checksums and the
   exact certification JSON; compare both extraction modes under the recorded
   event cap and each available version lock.
2. Save and compare event keys and selected float64 mass vectors; inspect the
   first differing element, if any.
3. If mass vectors match, compare histogram edges, bin membership, window
   masks, fit-eligibility/clipping and baseline arrays.
4. If baseline arrays match, compare residual arrays, locked template dot
   products, and the paired blocked-CV fold records and tie breaks.
5. Add a regression for the **confirmed** first cause; only then call the
   causal part of issue #9 complete.

No historical output is overwritten to force equality. The current small
numerical differences do not change the observed positive sign, pair
winner, or scientific classification: **reproducible pipeline-dependent
residual, Category C significance, physical origin unresolved**.
Issue #9 remains open until causal trace or a fully evidenced bounded
environment explanation is completed.
