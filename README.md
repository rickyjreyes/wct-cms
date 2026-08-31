# CMS WCT Validation

Blind/out-of-sample CMS NanoAOD cross-validation of a pre-specified residual signature previously identified outside CMS.

The repository currently implements the **opposite-sign dimuon invariant-mass channel**. Its core purpose is not to discover an arbitrary CMS oscillation; it is to ask whether a frequency/signature fixed from an external analysis survives in CMS under an independent detector and reconstruction chain.

## Tested signature

For residuals of a smooth dimuon mass spectrum,

\[
r(m)=\frac{N(m)-B(m)}{\sqrt{B(m)}},
\]

we test

\[
r(m)=c+a\cos[\omega\ln(m/m_0)]+b\sin[\omega\ln(m/m_0)].
\]

Equivalently,

\[
r(m)=c+A\cos[\omega\ln(m/m_0)-\phi].
\]

`--frozen-omega` is the primary replication test. The unrestricted omega scan is exploratory and gets a permutation-based global p-value.

## Current recorded replication result

The candidate frequency

```text
omega_m = 7.025825825825827
```

was identified in one Run2016H file and frozen before inspection of a second Run2016H file and before inspection of the preregistered Run2016G cross-period sample.

Both frozen replications pass the implemented tests:

| sample | amplitude | phase (rad) | Delta chi-square | local p (2 dof) | local one-sided Z |
|---|---:|---:|---:|---:|---:|
| Run2016H file-2 frozen replication | `0.9367120932` | `-0.3059910509` | `118.9148340` | `1.5065e-26` | `10.5990 sigma` |
| Run2016G cross-period frozen replication | `0.9348796605` | `-0.1567923197` | `115.8921026` | `6.8289e-26` | `10.4567 sigma` |

The Run2016G amplitude is within about `0.20%` of the Run2016H replication amplitude. Its phase is only about `1.84 degrees` from the original Run2016H discovery-file phase under the same `m0 = 1 GeV` convention.

For Run2016G, both the frozen residual-permutation test and the frozen parametric refit-bootstrap produced zero exceedances and hit their finite Monte Carlo floors:

```text
frozen permutation p = 1/1001 = 0.000999000999000999
frozen refit-bootstrap p = 1/501 = 0.001996007984031936
```

The approximately `10.5--10.6 sigma` values are **local analytic fixed-frequency diagnostics only**. They are not empirical/global significances and are conditional on the implemented background/residual model. The Run2016G residual field remains broader than an ideal Pearson field (`RMS = 1.7240`, `max |r| = 8.4725`), so model-adequacy and systematics tests remain necessary before any discovery-grade physical interpretation.

See:

- `docs/CMS_REPLICATION_FILE2_RESULT_2026-08-28.md` for the Run2016H frozen replication;
- `docs/CMS_RUN2016G_REPLICATION_FREEZE_2026-08-28.md` for the preregistered cross-period protocol;
- `docs/CMS_RUN2016G_RESULT_2026-08-31.md` for the Run2016G result and cross-sample comparison.

## Repository layout

```text
cms-wct-validation/
├── src/cms_wct/
│   ├── analysis.py       end-to-end pipeline
│   ├── background.py     histogram + robust smooth background
│   ├── cmsio.py          NanoAOD input + dimuon reconstruction
│   ├── signature.py      fixed-frequency and scanned statistics
│   ├── plots.py          diagnostic figures
│   ├── models.py         result dataclasses
│   └── cli.py            command-line interface
├── tests/                synthetic/unit tests
├── configs/              frozen-analysis template
├── data/                 input manifest template
├── docs/                 analysis protocol
├── .github/workflows/    CI
├── legacy_single_script.py
├── pyproject.toml
└── requirements.txt
```

## Install

```bash
python -m venv .venv
# Windows Git Bash:
source .venv/Scripts/activate
# Linux/macOS:
# source .venv/bin/activate

pip install -e .[dev]
pytest -q
```

## Input

Create `data/files.txt` containing one NanoAOD ROOT file or XRootD URL per line.

```text
root://.../file1.root
root://.../file2.root
```

ROOT files are intentionally ignored by git.

## Blind run

Freeze the externally predicted omega **before inspecting CMS output**:

```bash
cms-wct \
  --input data/files.txt \
  --output-dir results/dimuon_blind \
  --mass-min 2 \
  --mass-max 120 \
  --bins 350 \
  --log-bins \
  --muon-pt-min 4 \
  --muon-eta-max 2.4 \
  --tight-id \
  --fit-degree 7 \
  --omega-min 0.5 \
  --omega-max 80 \
  --omega-steps 3000 \
  --frozen-omega YOUR_PREEXISTING_VALUE \
  --permutations 2000 \
  --seed 20260827
```

For a quick pipeline test:

```bash
cms-wct --input data/files.txt --max-events 100000 --permutations 20 --frozen-omega YOUR_PREEXISTING_VALUE
```

## Outputs

Each run writes:

- `summary.json`
- `spectrum.csv`
- `omega_scan.csv`
- `permutation_global_max.csv`
- `permutation_frozen.csv` when a frozen omega is supplied
- `spectrum_background.png`
- `residuals.png`
- `omega_scan.png`
- `best_log_periodic_fit.png`
- `global_null_distribution.png`
- exact run configuration and input manifest

## Interpretation

A low p-value at the **frozen** frequency is the relevant replication statistic. A low global scan p-value means CMS contains some unusually strong frequency under this model, but it does not by itself replicate the external WCT prediction.

Before treating any result as physical evidence, run the systematic controls listed in `docs/ANALYSIS_PROTOCOL.md`, especially background-model stability, data-era splits, CMS simulation, and a second CMS channel.

## Current scope

This initial repository is deliberately focused on dimuons because the observable is clean, physically interpretable, and suitable for residual-spectrum cross-validation. Photon and jet channels should be added as independent modules rather than mixed into the first blind test.
