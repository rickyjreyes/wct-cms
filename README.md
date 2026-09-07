# CMS Run-2 Log-Periodic Dimuon Residual

Reproducible open-data analysis of a fixed-frequency residual in the CMS Run-2 opposite-sign dimuon invariant-mass spectrum, including frozen independent-file replication, cross-period replication, a prospective phase-locked holdout, and an independent adversarial statistical audit.

**Paper:** [Log-Periodic Dimuon Residual in CMS Open Data: Cross-Period Replication and a Prospective Phase-Locked Holdout](https://zenodo.org/records/22257067)  
**Independent audit:** [`docs/CMS_INDEPENDENT_ADVERSARIAL_AUDIT.md`](docs/CMS_INDEPENDENT_ADVERSARIAL_AUDIT.md)  
**Research:** [rickyjreyes.github.io](https://rickyjreyes.github.io)

> **Current status — 2026-09-07:** **C. Suggestive only / insufficiently robust.** The fixed-frequency residual reproduces, including on the prospectively phase-locked G2 file, but the previously quoted significance is not robust to reasonable smooth-background uncertainty. Alternative smooth generating means with no explicitly injected sinusoid exceed the observed H2/G1 pair score in **359–442 of 1,000 trials** when the original background-selection, fitting, and scoring pipeline is rerun. This does **not** prove that the physical spectrum contains no periodic component; it does invalidate a model-robust discovery-significance interpretation from the current pipeline alone.

## Independent adversarial audit

The audit independently reproduces the historical result before attacking its assumptions.

| Diagnostic | Audit result |
|---|---:|
| Reproduced H2/G1 selected pair score | `108.43746` |
| Reproduced G2 locked score | `126.27500` |
| Historical selected-background null | `0 / 10,000` exceedances |
| Spline `s=1` smooth generator | `442 / 1,000` exceedances |
| Poisson log-polynomial degree 12 | `359 / 1,000` exceedances |
| Penalized spline (16 knots, penalty 0.1) | `405 / 1,000` exceedances |
| Permissive degree-12 broad-frequency diagnostic | `918 / 1,000` exceedances |

The central obstruction is **background/signal identifiability**. Several defensible backgrounds predict held-out spectra better while removing the locked excess, but flexible backgrounds also absorb real injected waveforms. Therefore the audit stops short of calling the residual an artifact.

The historical analytic p-values and the `0/10,000` historical-generator result below remain useful **conditional diagnostics**, but they are not validated physical discovery probabilities.

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

The H2 and G1 amplitudes differ by only about `0.20%` under the historical background convention.

### G2 prospective phase-locked holdout

Before inspection of the G2 target file, the file-selection rule, file identity, frequency, phase, positive amplitude sign, event selection, mass range, binning, resonance masks, background model, null sizes, and random seed were frozen.

Observed historical result:

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
p_{\mathrm{refit}} = \frac{1}{501} \approx 1.9960\times10^{-3}.
$$

The extremely small analytic probability is a **fixed-waveform diagnostic conditional on the model**. The independent audit additionally shows that its physical significance is highly background-model dependent.

---

## Historical flexible-background kill test

The original replication sequence used a degree-7 Chebyshev continuum. A central concern was therefore whether background fitting or detrending could manufacture the frozen waveform, or whether a sufficiently flexible continuum could absorb it.

The historical pipeline implements a WCT-blind predictive-background selector over:

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

### Historical end-to-end smooth-null calibration

Each pseudoexperiment reruns:

1. smooth-background generation;
2. background-family selection;
3. continuum refitting;
4. residual construction;
5. the frozen waveform test.

Observed historical calibration:

$$
N_{\mathrm{exceed}} = 0/10{,}000
$$

with add-one Monte Carlo probability

$$
p_{\mathrm{MC}} = \frac{0+1}{10{,}000+1} = 9.9990\times10^{-5}.
$$

This shows that the observed score is unusual **conditional on that selected generating mean**. The independent audit demonstrates that this conclusion does not survive other reasonable smooth generating backgrounds: the same unchanged analysis pipeline produces exceedance fractions of roughly `0.36–0.44` under three alternatives.

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

relative to the historical smooth-null 95th-percentile threshold.

The audit confirms that flexible backgrounds can also absorb injected signal, especially at higher flexibility. This is why low $Q$ under a flexible background is not, by itself, proof of an artifact.

Run the canonical historical background-kill pipeline with:

```bash
python scripts/run_cms_background_kill.py \
  --null-trials 10000 \
  --injection-trials 1000 \
  --injection-amplitudes 0.25 0.5 0.75 1.0
```

Run the independent adversarial suite with:

```bash
python -m pip install -r experiments/cms_audit/requirements.txt
python -m pip install -e '.[dev]'
python experiments/cms_audit/fetch_inputs.py
OPENBLAS_NUM_THREADS=1 python experiments/cms_audit/reproduce_initial.py
OPENBLAS_NUM_THREADS=1 python experiments/cms_audit/run_spectrum_attacks.py
OPENBLAS_NUM_THREADS=1 python experiments/cms_audit/run_failure_certificate.py --out results/cms_audit_certificate_clean
OPENBLAS_NUM_THREADS=1 python experiments/cms_audit/run_additional_controls.py
OPENBLAS_NUM_THREADS=1 python experiments/cms_audit/run_model_averaging.py
OPENBLAS_NUM_THREADS=1 python experiments/cms_audit/run_event_controls.py
python -m pytest -q
```

---

## What the current result establishes

The combined historical analysis and independent audit support the following narrower empirical statements:

1. an interior log-frequency selected in one certified Run2016H file reproduced at the frozen frequency in an independent Run2016H file;
2. the same frozen frequency reproduced in a separately preregistered Run2016G cross-period test;
3. a previously unused Run2016G file supported the already-frozen frequency, phase, and positive amplitude sign under the historical background procedure;
4. the fixed-frequency residual is not a one-bin, simple binning, or basic numerical-arithmetic accident;
5. the historical selected-background null produces `0/10,000` exceedances, but alternative defensible smooth generating means produce `359–442/1,000` exceedances through the unchanged pipeline;
6. flexible backgrounds can absorb a genuine injected waveform, so the audit does not establish that the physical residual is an artifact.

The correct current interpretation is **reproducible fixed-frequency structure under the historical analysis assumptions, with unresolved background/signal identifiability**.

---

## What it does not establish

The current result does not by itself show that:

- WCT is the unique physical cause;
- the residual is a new particle or resonance;
- CMS detector or reconstruction effects cannot generate it;
- trigger/selection/acceptance structure cannot generate it;
- correlated detector systematics are negligible;
- Standard Model continuum, resonance tails, or interference cannot generate it;
- the analytic local tail is a calibrated physical p-value;
- the historical `0/10,000` result is robust across reasonable background uncertainty;
- the empirical tail probability is $>5\sigma$;
- the CMS result and results in other physical domains are statistically independent evidence for one common mechanism.

There is currently **no unique validated global p-value** for the CMS claim in this repository.

---

## Next falsification priorities

The single highest-value test is an independently constrained prospective holdout:

1. **Preregister an unseen CMS subset with an independently constrained background** — use an efficiency-matched control sample or validated detector-folded Standard Model prediction; freeze frequency, phase, trigger plateau, masks, nuisance/background procedure, and injection-recovery acceptance before opening the target.
2. **Trigger and reconstruction efficiency controls** — test whether known efficiency structure projects onto the frozen waveform.
3. **Correlated detector/systematic nulls** — replace independent smooth-Poisson pseudoexperiments with justified correlated uncertainty models.
4. **Standard Model and resonance/interference controls** — propagate broad continuum, resonance tails, and interference models through the identical residual pipeline.
5. **Acceptance and selection tests** — stress muon kinematics, IDs, masks, run subdivisions, and detector-era structure.
6. **Independent detector replication** — test the frozen observable/signature with ATLAS or another genuinely independent detector chain where compatible data exist.

More trials around the same historical fitted mean do **not** resolve the background/signal identifiability problem.

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
│   ├── background_kill.py      historical end-to-end null + injection tests
│   ├── cmsio.py                NanoAOD input + dimuon reconstruction
│   ├── signature.py            fixed-frequency and scanned statistics
│   ├── locked.py               fixed-frequency/fixed-phase directional tests
│   ├── significance.py         Monte Carlo resolution and exact tail bounds
│   ├── plots.py                diagnostic figures
│   ├── models.py               result dataclasses
│   └── cli.py                  command-line interface
├── experiments/cms_audit/      independent adversarial audit implementation
├── audit/2026-09-06/           compact archived audit evidence and summaries
├── scripts/
├── tests/
├── configs/
├── data/
├── docs/
├── .github/workflows/
├── legacy_single_script.py
├── pyproject.toml
└── requirements.txt
```

ROOT inputs and regenerable per-chunk Monte Carlo audit checkpoints are intentionally ignored by git. The audit keeps compact manifests, aggregate summaries, regression witnesses, tables, and plots under version control.

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

The sharpest historical holdout script freezes frequency, phase, and positive amplitude sign:

```text
scripts/run_phase_locked_period.py
```

Canonical G2 result record:

```text
docs/CMS_RUN2016G_FILE2_PHASE_LOCK_RESULT_2026-08-31.json
```

---

## Background kill test

Quick/default historical diagnostic:

```bash
python scripts/run_cms_background_kill.py
```

Deep run matching the historical selected-background calibration:

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

Zero exceedances in 10,000 trials are nowhere near enough to resolve this tail directly. More importantly, the independent audit shows that increasing trial count around the same historical fitted mean would not solve the larger background-model uncertainty.

The repository retains the historical direct-Monte-Carlo planning document in:

```text
docs/EMPIRICAL_5SIGMA_PROTOCOL_2026-08-31.md
```

For zero exceedances:

| criterion | required trials |
|---|---:|
| add-one numerical floor reaches $5\sigma$ p scale | `3,488,555` |
| exact one-sided 95% upper bound reaches threshold | `10,450,778` |
| exact one-sided 99% upper bound reaches threshold | `16,065,391` |

Those trial counts matter only after the composite background/nuisance model is independently justified. No combined H/G/G2 sigma is reported by multiplying p-values or adding Z values.

---

## Interpretation hierarchy

Keep three claims separate:

1. **Empirical:** the declared historical pipeline reproducibly returns a fixed-frequency residual across H2, G1, and the phase-locked G2 holdout.
2. **Statistical robustness:** unresolved; reasonable smooth-background alternatives remove the positive component and generate historical-sized scores routinely in the unchanged pipeline.
3. **Physical attribution:** whether the recurring structure is signal, detector/acceptance background, Standard Model structure, or another mechanism remains open.

The repository now preserves both the original positive replication chain and the independent failure of its model-robust significance interpretation.