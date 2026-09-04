# CMS Run-2 Log-Periodic Dimuon Residual

Reproducible open-data analysis of a log-periodic residual in the CMS Run-2 opposite-sign dimuon invariant-mass spectrum, with frozen independent-file replication, cross-period replication, a prospective phase-locked holdout, and an end-to-end flexible-background kill test.

**Paper:** [Log-Periodic Dimuon Residual in CMS Open Data: Cross-Period Replication and a Prospective Phase-Locked Holdout](https://zenodo.org/records/22257067)  
**Research:** [rickyjreyes.github.io](https://rickyjreyes.github.io)

> **Current status:** positive replication under the implemented CMS pipeline, including a phase-locked holdout and a 10,000-trial end-to-end smooth-background adversarial calibration. This is **not yet a discovery-grade physical anomaly** and does **not uniquely establish Wave Confinement Theory (WCT)**. Detector, reconstruction, acceptance, correlated-systematic, Standard Model, resonance/interference, and independent-detector controls remain open.

---

## What this repository tests

The primary observable is the inclusive opposite-sign dimuon invariant mass

$$
m_{\mu\mu}
$$

analyzed in the logarithmic coordinate

$$
x = \ln\left(\frac{m_{\mu\mu}}{1\,\mathrm{GeV}}\right).
$$

After fitting a smooth continuum background $B(m)$, Pearson-like residuals are defined as

$$
r(m)=\frac{N(m)-B(m)}{\sqrt{B(m)}}.
$$

The tested residual model is

$$
r(m)=c+a\cos(\omega x)+b\sin(\omega x)
$$

or equivalently

$$
r(m)=c+A\cos(\omega x-\phi).
$$

The frozen CMS frequency is

$$
\omega_{\mathrm{CMS}} = 7.025825825825827
$$

in $\ln(m_{\mu\mu}/1\,\mathrm{GeV})$.

The important distinction is chronological:

- WCT motivated a **pre-existing prediction class** of log-periodic collider structure before this CMS analysis;
- the specific numerical CMS frequency $\omega_{\mathrm{CMS}} = 7.025825825825827$ was selected in the first certified Run2016H discovery file;
- that numerical value was then frozen before the subsequent independent-file and cross-period tests.

Do not relabel this CMS frequency as $k \sim 9.7$. That value belongs to a different observable/coordinate in the GWTC analysis.

---

## Evidence chain

The analysis progressively removes fitting freedom.

| Stage | Dataset | What was free? | Amplitude | Phase (rad) | $\Delta\chi^2$ |
|---|---|---|---:|---:|---:|
| **H1 discovery** | Run2016H file 1 | frequency + phase | `0.7543` | `-0.1890` | `75.76` |
| **H2 frozen replication** | independent Run2016H file 2 | phase only at frozen frequency | `0.9367121` | `-0.3059911` | `118.9148` |
| **G1 cross-period replication** | preregistered Run2016G file 1 | phase only at frozen frequency | `0.9348797` | `-0.1567923` | `115.8921` |
| **G2 phase-locked holdout** | previously unused Run2016G file 2 | positive amplitude only; frequency + phase frozen | `0.9708618` | `-0.2313917` frozen | `126.2832` |

The H2 and G1 amplitudes differ by only about `0.20%`.

### G2 prospective phase-locked holdout

Before inspection of the G2 target file, the file-selection rule, file identity, frequency, phase, positive amplitude sign, event selection, mass range, binning, resonance masks, background model, null sizes, and random seed were frozen.

Observed result:

$$
A = 0.9708617746
$$

$$
\Delta\chi^2 = 126.2832399542
$$

$$
p_{\mathrm{analytic}} = 1.3329276765\times 10^{-29}
$$

for the one-sided fixed-waveform analytic diagnostic.

Finite Monte Carlo ensembles gave zero exceedances:

```text
residual permutations:                 0 / 1000
end-to-end Poisson background refits: 0 / 500
```

Therefore the empirical probabilities from those ensembles are limited by their Monte Carlo floors:

$$
p_{\mathrm{perm}} = \frac{1}{1001} \approx 9.9900\times10^{-4}
$$

$$
p_{\mathrm{refit}} = \frac{1}{501} \approx 1.9960\times10^{-3}
$$

The extremely small analytic probability is a **fixed-waveform diagnostic conditional on the model**. It is not an empirical $>5\sigma$ claim.

---

## Flexible-background kill test — completed

The original replication sequence used a degree-7 Chebyshev continuum. A central concern was therefore whether background fitting or detrending could manufacture the frozen waveform, or whether a sufficiently flexible continuum could absorb it.

The repository implements a WCT-blind predictive-background selector over:

```text
Chebyshev degrees 5..12
Bernstein degrees 5, 7, 9, 12
smoothing splines with factors 0.5, 1, 2
```

Backgrounds are ranked by blocked held-out Poisson deviance **without using the WCT frequency, phase, amplitude, or test statistic**.

The signal-independent procedure selected:

```text
spline_s2
```

The conservative H2-G1 pair statistic is

$$
T_{\mathrm{pair}}
=\min\!\left(\Delta\chi^2_{\mathrm{locked,H2}},\Delta\chi^2_{\mathrm{locked,G1}}\right)
=108.4978.
$$

### End-to-end smooth-null calibration

Each pseudoexperiment reruns:

1. smooth-background generation;
2. background-family selection;
3. continuum refitting;
4. residual construction;
5. the frozen waveform test.

Observed calibration:

$$
N_{\mathrm{exceed}} = 0/10{,}000
$$

with add-one Monte Carlo probability

$$
p_{\mathrm{MC}} = \frac{0+1}{10{,}000+1} = 9.9990\times10^{-5}.
$$

This materially reduces the specific explanation that the H/G structure is produced solely by the tested smooth-continuum selection or detrending procedure.

It does **not** calibrate all detector or Standard Model systematics.

### Injection / recovery

The same end-to-end pipeline was tested after injecting the frozen waveform at amplitudes

$$
A_{\mathrm{inj}} \in \{0.25,\,0.50,\,0.75,\,1.00\}.
$$

Across that range, the selected flexible-background pipeline retained approximately

$$
R_A \approx 0.76\text{--}0.81
$$

of the injected waveform amplitude and achieved approximately

$$
\mathrm{power} \approx 0.93\text{--}0.98
$$

relative to the smooth-null primary-score threshold.

This matters because a flexible background that simply absorbs every injected waveform would not provide an informative falsification test.

Run the canonical adversarial pipeline with:

```bash
python scripts/run_cms_background_kill.py \
  --null-trials 10000 \
  --injection-trials 1000 \
  --injection-amplitudes 0.25 0.5 0.75 1.0
```

Reusable implementation:

```text
src/cms_wct/background_kill.py
```

Regression tests:

```text
tests/test_background_kill.py
```

---

## What the current result establishes

Under the implemented CMS pipeline, the data support the following narrower empirical statements:

1. an interior log-frequency selected in one certified Run2016H file reproduced at the frozen frequency in an independent Run2016H file;
2. the same frozen frequency reproduced in a separately preregistered Run2016G cross-period test;
3. a previously unused Run2016G file supported the already-frozen frequency, phase, and positive amplitude sign;
4. a signal-independent flexible-background selector retained a large conservative H2-G1 locked statistic;
5. zero of 10,000 complete smooth-Poisson null pipelines produced an equally large pair statistic;
6. injection/recovery shows that the selected flexible continuum retains substantial sensitivity to the frozen waveform.

These statements support **reproducible log-periodic structure in the analyzed CMS observable under the declared pipeline and null models**.

---

## What it does not establish

The current result does not by itself show that:

- WCT is the unique physical cause;
- the residual is a new particle or resonance;
- CMS detector or reconstruction effects cannot generate it;
- trigger/selection/acceptance structure cannot generate it;
- correlated detector systematics are negligible;
- Standard Model continuum, resonance tails, or interference cannot generate it;
- the empirical tail probability is $>5\sigma$;
- the CMS result and results in other physical domains are statistically independent evidence for one common mechanism.

The remaining conventional explanation is a reproducible CMS/Standard-Model/acceptance/reconstruction structure not represented by the current smooth-Poisson null ensemble.

---

## Next falsification priorities

The highest-value next tests are physical/systematic controls, not another same-model sigma calculation:

1. **Trigger and reconstruction efficiency controls** — test whether known efficiency structure projects onto the frozen waveform.
2. **Correlated detector/systematic nulls** — replace independent smooth-Poisson pseudoexperiments with justified correlated uncertainty models.
3. **Standard Model and resonance/interference controls** — propagate broad continuum, resonance tails, and interference models through the identical residual pipeline.
4. **Acceptance and selection tests** — stress muon kinematics, IDs, masks, run subdivisions, and detector-era structure.
5. **Independent detector replication** — test the frozen observable/signature with ATLAS or another genuinely independent detector chain where compatible data exist.
6. **Independent CMS channels** — add clean channels as separate modules rather than mixing them into the dimuon discovery pipeline.

The ranked adversarial program is documented in:

```text
docs/CMS_KILL_TEST_PROGRAM_2026-09-01.md
```

The background-kill implementation is documented in:

```text
docs/CMS_BACKGROUND_KILL_TEST_IMPLEMENTATION_2026-09-01.md
```

---

## Frequency conventions

This repository uses

$$
x_{\mathrm{CMS}} = \ln\left(\frac{m_{\mu\mu}}{1\,\mathrm{GeV}}\right),
\qquad
\omega_{\mathrm{CMS}} = 7.025825825825827.
$$

The mapped LHCb request-48 work in `rickyjreyes/LHC` uses $\ln(q^2)$. Since $q^2=m^2$,

$$
k_{\mathrm{LHCb}} = 3.512912912912913,
\qquad
\omega_{\mathrm{CMS}} = 2k_{\mathrm{LHCb}}.
$$

Raw numerical frequencies from different logarithmic coordinates must not be compared without the coordinate conversion.

---

## Repository layout

```text
wct-cms/
├── src/cms_wct/
│   ├── analysis.py             end-to-end base pipeline
│   ├── background.py           base smooth background
│   ├── background_families.py  Chebyshev/Bernstein/spline fits
│   ├── background_cv.py        WCT-blind blocked predictive selection
│   ├── background_kill.py      end-to-end adversarial null + injection tests
│   ├── cmsio.py                NanoAOD input + dimuon reconstruction
│   ├── signature.py            fixed-frequency and scanned statistics
│   ├── locked.py               fixed-frequency/fixed-phase directional tests
│   ├── significance.py         Monte Carlo resolution and exact tail bounds
│   ├── plots.py                diagnostic figures
│   ├── models.py               result dataclasses
│   └── cli.py                  command-line interface
├── scripts/
│   ├── run_phase_locked_period.py
│   ├── run_hg_background_cv.py
│   ├── run_cms_background_kill.py
│   └── plan_empirical_5sigma.py
├── tests/
│   └── test_background_kill.py
├── configs/
├── data/
├── docs/
├── .github/workflows/
├── legacy_single_script.py
├── pyproject.toml
└── requirements.txt
```

ROOT inputs are intentionally ignored by git.

---

## Install

```bash
python -m venv .venv

# Windows Git Bash
source .venv/Scripts/activate

# Linux/macOS
# source .venv/bin/activate

pip install -e .[dev]
pytest -q
```

---

## Input

Create `data/files.txt` containing one NanoAOD ROOT file or XRootD URL per line:

```text
root://.../file1.root
root://.../file2.root
```

---

## Base frozen-frequency run

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
  --frozen-omega 7.025825825825827 \
  --permutations 2000 \
  --seed 20260827
```

The unrestricted omega scan is exploratory. The scientific replication statistic is the statistic evaluated at the frequency frozen before the target sample was inspected.

---

## Phase-locked prospective test

The sharpest holdout script freezes frequency, phase, and positive amplitude sign:

```text
scripts/run_phase_locked_period.py
```

Canonical G2 result record:

```text
docs/CMS_RUN2016G_FILE2_PHASE_LOCK_RESULT_2026-08-31.json
```

---

## Background kill test

Quick/default diagnostic:

```bash
python scripts/run_cms_background_kill.py
```

Deep run matching the current adversarial calibration:

```bash
python scripts/run_cms_background_kill.py \
  --null-trials 10000 \
  --injection-trials 1000 \
  --injection-amplitudes 0.25 0.5 0.75 1.0
```

Default outputs are written under:

```text
results/cms_background_kill/
```

including:

```text
selection_freeze.json
cv_scores.csv
spurious_null_trials.csv
spurious_null_summary.json
absorption_matrix.csv
injection_trials.csv
injection_summary.json
summary.json
```

---

## Empirical $>5\sigma$ protocol

For a one-sided Gaussian convention,

$$
5\sigma \iff p = 2.866515718791946\times10^{-7}.
$$

Zero exceedances in 10,000 trials are nowhere near enough to resolve this tail directly.

The repository therefore keeps a separate discovery-grade direct-Monte-Carlo protocol in:

```text
docs/EMPIRICAL_5SIGMA_PROTOCOL_2026-08-31.md
```

For zero exceedances:

| criterion | required trials |
|---|---:|
| add-one numerical floor reaches $5\sigma$ p scale | `3,488,555` |
| exact one-sided 95% upper bound reaches threshold | `10,450,778` |
| exact one-sided 99% upper bound reaches threshold | `16,065,391` |

The preferred direct-Monte-Carlo gate is the **95% exact upper-bound criterion**, together with the predeclared systematic-control envelope.

```bash
python scripts/plan_empirical_5sigma.py
```

No combined H/G/G2 sigma is reported by multiplying p-values or adding Z values. A combined significance requires a joint statistic and joint null procedure frozen in advance.

---

## Interpretation hierarchy

Keep three claims separate:

1. **Empirical:** a reproducible log-periodic residual is present in the analyzed CMS dimuon observable under the implemented pipeline.
2. **Theoretical consistency:** the residual is consistent with the pre-existing WCT prediction class of log-periodic collider structure.
3. **Physical attribution:** whether the residual is caused by WCT rather than detector, Standard Model, statistical, or other discrete-scale mechanisms remains open.

The repository is designed to make claim 1 increasingly difficult to explain as an analysis artifact while keeping claim 3 explicitly falsifiable.
