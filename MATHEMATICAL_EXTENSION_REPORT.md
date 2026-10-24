# Mathematical extension: applicability audit

## Scope

This extension tests which mathematical tools are identifiable from the
bundled single-section spatial snapshot. It is deliberately an audit, not a
claim that a prestigious mathematical object explains tumour biology.

The reproducible implementation is
[`src/mathematical_applicability_audit.py`](./src/mathematical_applicability_audit.py).
It builds a symmetric six-nearest-neighbour graph for each cached cohort,
computes low-frequency normalized-Laplacian eigenvalues, a phenotype
Dirichlet energy, and a **truncated spectral zeta**

`zeta_L^(k)(s) = sum(lambda_i^(-s))`

over the retained positive eigenvalues. This is not the Riemann zeta function
and not a full graph zeta invariant. It is a descriptive graph-spectrum
summary whose value depends on the graph construction and truncation.

## Results

The machine-readable output is
[`results/MATHEMATICAL_APPLICABILITY_AUDIT.json`](./results/MATHEMATICAL_APPLICABILITY_AUDIT.json).

| Cohort | Nodes | Edges | Components | Truncated zeta s=1 | Phenotype Moran-like |
|---|---:|---:|---:|---:|---:|
| Breast | 4,992 | 15,227 | 1 | 7,841.25 | 0.709 |
| Melanoma | 3,458 | 10,785 | 2 | 5,860.84 | 0.729 |
| Lung | 3,858 | 11,859 | 1 | 6,704.51 | 0.731 |

These values are **diagnostics only**. The melanoma graph is disconnected,
and the graph spectrum changes with `k`, boundary handling, and coordinate
quality. The positive association between phenotype signal and neighbouring
spots is compatible with spatial clustering; it does not identify a game
interaction, a payoff, a clone, or a causal mechanism.

## What can and cannot be applied

### Identifiable with current data

- pair-correlation and Ripley-type summaries for spatial clustering;
- graph Laplacian and low-frequency spectral diagnostics;
- patient-level bootstrap and spatially blocked validation;
- conditional-information analyses, only with bias-corrected estimators and
  spatially structured permutations.

### Requires new measurements

- Z-transform, Koopman operators, SDEs, or transition-rate estimation:
  repeated measurements or a discrete time series;
- Pontryagin control:
  intervention response, a calibrated state transition model, toxicity and
  treatment constraints;
- genetic selection:
  mutation/CNV/lineage measurements or a matched single-cell reference.

### Not justified as a current biological method

- Riemann zeta: the snapshot supplies no observable, likelihood, or null
  hypothesis that would make `zeta(s)` identifiable;
- quantum, relativistic, or Millennium Prize analogies: they add no
  distinguishable prediction at this measurement scale.

## Interpretation boundary

The current audit supports a statement about spatial graph structure only.
It does not overturn the existing patient-held-out result: the game
advantage remains small and its bootstrap interval includes zero. A future
claim would require a pre-specified statistic, an appropriate spatial null,
patient-held-out improvement over composition and spatial baselines, and
stability across marker panels and graph definitions.

The literature review used for method selection reports Ripley/pair-correlation
and graph-based spatial analyses as established spatial-transcriptomics tools;
zeta-based extensions remain niche and descriptive in this setting.
