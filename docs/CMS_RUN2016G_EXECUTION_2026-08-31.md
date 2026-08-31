# CMS Run2016G frozen replication execution

This note records the executable command for the preregistered Run2016G cross-period replication defined in `CMS_RUN2016G_REPLICATION_FREEZE_2026-08-28.md`.

No analysis choice is changed here. The command below maps the already-frozen protocol directly onto the current `cms-wct` CLI.

## Input files

- Input manifest: `data/files_run2016g.txt`
- First-file sample: `data/cms_run2016g/30522/05DD095C-F6C3-9A4F-9FB3-348A5A6403D5.root`
- Golden JSON: `data/cert/14221/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON_MuonPhys.txt`
- Output directory: `results/run2016g_frozen`

## Canonical execution

The repository root now contains `run_run2016g_frozen.sh`, which executes:

```bash
cms-wct \
  --input data/files_run2016g.txt \
  --golden-json data/cert/14221/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON_MuonPhys.txt \
  --output-dir results/run2016g_frozen \
  --mass-min 2 \
  --mass-max 120 \
  --bins 350 \
  --log-bins \
  --muon-pt-min 4 \
  --muon-eta-max 2.4 \
  --tight-id \
  --max-events 100000 \
  --fit-degree 7 \
  --mask-window 2.9:3.3 \
  --mask-window 3.55:3.85 \
  --mask-window 8.5:11.5 \
  --mask-window 80:100 \
  --omega-min 3.1 \
  --omega-max 80 \
  --omega-steps 1000 \
  --frozen-omega 7.025825825825827 \
  --permutations 1000 \
  --parametric-bootstrap 500 \
  --seed 20260827
```

The repository default values for fit iterations (`6`) and clip sigma (`3.5`) are deliberately not overridden because the frozen protocol specifies the repository defaults.

## Interpretation lock

The primary result is the fixed-frequency result at exactly `omega_m = 7.025825825825827`. The unrestricted best-frequency scan remains secondary. Do not alter the frequency, mass windows, event cap, binning, cuts, background degree, or null counts after seeing the Run2016G result and still call the altered run the preregistered replication.

Record at minimum from `results/run2016g_frozen/summary.json`:

- `frozen_scan.amplitude`
- `frozen_scan.phase`
- `frozen_scan.delta_chi2`
- `frozen_permutation_p`
- `frozen_parametric_bootstrap_p`
- `residual_rms`
- `residual_max_abs`
- `best_scan_at_boundary`
- `best_cycles_across_span`
