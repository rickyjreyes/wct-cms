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

The preregistered Run2016H file-2 replication at the frozen frequency

```text
omega_m = 7.025825825825827
```

returned

```text
delta chi-square = 118.91483403371683
local chi-square p (2 dof) = 1.506509523217018e-26
one-sided Gaussian-equivalent Z = 10.5989632046 sigma
```

Thus the frozen-frequency analytic local diagnostic is approximately **10.6 sigma**.

This number is deliberately labeled **local and diagnostic only**. It is the Gaussian-equivalent mapping of the analytic fixed-frequency chi-square p-value; it is not the residual-permutation significance, the parametric-bootstrap significance, a look-elsewhere-corrected/global significance, or a claim of 10.6-sigma physical discovery. The empirical null tests in this finite run reached their Monte Carlo floors (`1/1001` for residual permutation and `1/501` for the refit bootstrap), and model-adequacy/systematic controls remain required.

See `docs/CMS_REPLICATION_FILE2_RESULT_2026-08-28.md` for the frozen result, phase comparison, empirical-null results, and model-adequacy warning.

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
