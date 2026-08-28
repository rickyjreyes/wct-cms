# Analysis protocol

## Primary question

Does a signature specified before inspection of CMS data appear in an independent CMS observable?

The primary model tested by this repository is

\[
r(m)=c+A\cos\left(\omega\ln(m/m_0)-\phi\right)+\epsilon,
\]

where `omega` should be frozen from the external WCT analysis before the CMS result is inspected.

## Primary vs exploratory result

**Primary replication test:** `--frozen-omega`.

The frequency is not fit to CMS. Only amplitude, phase, and intercept are fit. The repository reports the fixed-frequency improvement and a permutation p-value.

**Exploratory diagnostic:** omega scan.

The maximum scan statistic is corrected using the distribution of the maximum statistic obtained in residual permutations. It must not be reported as though it were the pre-registered replication test.

## Controls required before a physics claim

1. Repeat with multiple smooth-background families/degrees.
2. Vary binning without choosing the variation that maximizes significance.
3. Repeat across data-taking eras.
4. Run CMS simulation/background samples through the identical pipeline.
5. Test known detector/reconstruction resonances and masked windows.
6. Test at least one independent CMS channel, preferably photons after dimuons.
7. Apply the same dimensionless mapping used to compare GWTC, LHC, JUNO and photodiode data.
8. Record all failed as well as successful runs.

A CMS match by itself is not evidence that WCT caused the structure. The intended value is independent cross-validation under a frozen prediction and systematic controls.
