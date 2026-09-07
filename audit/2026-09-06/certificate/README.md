# CMS audit failure certificate

This directory intentionally keeps only compact, review-stable certificate outputs in Git:

- `manifest.json` — frozen source commit, observed threshold, analysis configuration, case counts, and deterministic seed rules.
- `summary.json` — aggregate exceedance counts, tail estimates, quantiles, amplitudes, and injection-recovery summaries.

The audit computation itself runs 18,000 pair-pipeline trials. Individual per-chunk checkpoint JSON files are regenerable working artifacts and are intentionally excluded from version control to keep the repository and pull-request history compact.

Regenerate the full certificate from the repository root with:

```bash
OPENBLAS_NUM_THREADS=1 python experiments/cms_audit/run_failure_certificate.py \
  --out results/cms_audit_certificate_clean
```

The deterministic seeds and case counts required to reproduce the aggregate certificate are recorded in `manifest.json`. The seeded smooth-background counterexample is additionally locked by `test_smooth_background_seeded_failure_witness` in `tests/test_cms_adversarial_audit.py`.
