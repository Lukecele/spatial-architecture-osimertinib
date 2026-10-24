# Rigorous bioinformatic protocol (October 2026)

## Primary question

For each tumour and platform separately, does a tumour-adapted spatial context
add reproducible information beyond the focal molecular state?

The protocol does not require the answer to be positive and does not require a
game-theoretic interpretation.

## Boolean filters

Every spot enters analysis through an auditable boolean mask:

```text
KEEP =
  tissue_present
  AND total_counts >= min_counts
  AND detected_features >= min_features
  AND total_counts <= max_counts (when pre-specified)
  AND detected_features <= max_features (when pre-specified)
  AND classification NOT IN {Artefact, Necrosis}
```

Thresholds must be frozen per platform/cohort before comparing outcomes. A
filter is not allowed to depend on the held-out outcome or on which model wins.
The implementation is [`src/rigorous_stats.py`](./src/rigorous_stats.py).

## Matrix layer

The matrix layer is explicit and versionable:

1. raw non-negative count matrix;
2. library-size normalised `log1p` CPM;
3. module or cell-composition matrix;
4. optional feature-wise z-score matrix;
5. focal and neighbourhood design matrices;
6. patient-level summary/effect matrix.

The z-transform is a numerical standardisation, not a biological law and not a
claim that a gene follows a special “z function”. Matrix orientation is always
recorded as genes-by-spots or spots-by-features.

## Statistical hierarchy

### Descriptive and spatial

- spatial graph and coordinate audit;
- Moran/Geary or variogram diagnostics;
- radius and kNN sensitivity;
- restricted spatial permutation nulls.

### Predictive

- focal-only baseline;
- focal plus context;
- smoothing and region-only controls;
- grouped patient-held-out evaluation;
- nested scale selection where scale is learned.

### Inference

- patient/donor is the biological replicate;
- bootstrap or hierarchical intervals are patient-level;
- Benjamini-Hochberg correction is used for families of hypotheses;
- spot-level p-values are never presented as patient-level evidence.

### Bayesian sensitivity

The lightweight normal-normal posterior in
[`src/rigorous_stats.py`](./src/rigorous_stats.py) can summarise patient-level
effects with an explicit prior. Prior sensitivity must be reported. PyMC or
NumPyro hierarchical models are optional upgrades only when the number of
patients, likelihood, and variance structure justify them.

## AI and advanced bioinformatics

Candidate tools are not automatically primary analyses:

- **cell2location, RCTD, stereoscope, Tangram**: only with an appropriate
  matched single-cell reference and uncertainty audit;
- **Squidpy, Giotto, Seurat spatial, BayesSpace, SpaGCN/STAGATE**: domain,
  neighbourhood and spatial-variable diagnostics, benchmarked against simple
  baselines;
- **Cellpose/StarDist**: only when image resolution and segmentation quality
  support cell-level claims;
- **scVI/scANVI**: integration or latent representations, never used to hide
  batch or patient leakage;
- **GNNs and foundation models such as Novae**: P2/P3 exploratory models
  requiring nested grouped validation, leakage checks and a sample-size
  justification;
- **ligand-receptor tools**: hypothesis generation, not causal communication;
- **PyMC/NumPyro**: hierarchical uncertainty when the data support it.

Recent method families such as SpaInGNN, STING, AugGCL and Novae are recorded
in the expanded inventory, but no high-capacity model is allowed to replace the
P0 baselines.
The source register is maintained in
[`CURRENT_METHODS_EVIDENCE_2026.md`](./CURRENT_METHODS_EVIDENCE_2026.md).

## Mandatory negative controls

- shuffled neighbourhood labels;
- sample-label permutation;
- radius/kNN alternative;
- region-only or coordinate-only control;
- QC-only control;
- marker-panel sensitivity;
- leave-one-patient-out replication.

## Decision rule

An effect enters the primary claim only if it:

1. beats a relevant simple baseline under grouped validation;
2. survives pre-specified sensitivity analyses;
3. has patient-level uncertainty that excludes the null or is clearly reported
   as inconclusive;
4. is not explained by QC, region, depth, or platform;
5. has an interpretation no stronger than the measured data.

## Current status

The current public breast and lung replications are negative/uncertain for the
fixed 150 µm model. This is not a reason to select a more complex AI model
post hoc. It is a reason to run P0 diagnostics and tumour-specific scale and
representation audits before any P2/P3 model.
