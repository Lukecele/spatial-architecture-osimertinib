# Interdisciplinary validation framework

## Purpose

This document extends the integrated InSituGame study without claiming a
literal "theory of everything". The object remains an observational,
spatially indexed model: local molecular state and neighborhood composition
are related to a proliferation proxy. Every added theory must produce a
testable model, a required data input, and a falsifiable output.

The original project and the public paired pilot remain separate evidence
layers. No new method below is considered validated merely because it has a
name from another discipline.

## Priority 0: establish that the measured interaction is real

| Area | Implementation | Required data | Falsifiable result | Main limitation |
|---|---|---|---|---|
| Measurement error and deconvolution | Compare at least two reference-based methods, such as cell2location, RCTD or Stereoscope; propagate posterior uncertainty into the payoff fit. | Matched single-cell/single-nucleus reference, raw spatial counts, platform and spot metadata. | Game parameters remain directionally stable across references and deconvolution methods; posterior predictive checks reproduce counts. | A transcriptomic state is not automatically a genetic clone or cell count. |
| Spatial statistics | Compare the game model with composition-only, density-only, distance-decay, neighborhood and spatial-error/random-effects models. Use spatially structured permutations and multiple-testing control. | Coordinates, boundaries, histology masks, adjacency graph and negative controls. | Interaction terms improve held-out spatial-block prediction and disappear under appropriate nulls. | Spatial autocorrelation can mimic interaction. |
| Hierarchical models | Use patient, specimen, section, region, batch and platform as separate levels; validate with patient-held-out splits. | Multiple independent patients and complete provenance metadata. | Estimate between-patient variance and out-of-patient performance instead of spot-only significance. | Many spots do not equal many biological replicates. |
| Uncertainty and sensitivity | Use block bootstrap or Bayesian hierarchical models; vary gene signatures, niche covariates and spatial graphs. | Replicate sections, measurement-error estimates and alternative marker sets. | Game class and equilibria remain stable under plausible perturbations with calibrated intervals. | Narrow intervals can still reflect model misspecification. |

## Priority 1: connect association to mechanism

| Area | Implementation | Required data | Falsifiable result | Main limitation |
|---|---|---|---|---|
| Causal inference | Define treatment/perturbation, outcome, confounders and estimand before analysis. Use randomized ex vivo perturbation where possible; otherwise target-trial emulation, g-computation or doubly robust estimation. Treat neighborhood interference explicitly. | Longitudinal or perturbation data, exposure timing, baseline covariates and outcomes. | Pre-specified treatment or mediation effect replicates in an independent experiment. | Cross-sectional adjustment does not establish causality. |
| Population genetics | Combine spatial expression with CNV, mutation, methylation or lineage-barcode data. Estimate clone frequencies and selection while testing whether states are reversible. | Spatial DNA/lineage data or validated genomic references. | Putative T_R and T_S states map to stable genetic/lineage groups, or are relabeled as phenotypic states. | Evolutionary-game coefficients are not automatically genetic fitness. |
| Dynamical systems and control | Compare replicator, Lotka-Volterra, ecological logistic and state-space models. Use model-predictive control only after transition rates are learned from time series. | Serial biopsies, treatment exposure, tumor burden and repeated state measurements. | Prospective prediction after a dose change; policy wins in a pre-registered experiment. | A spatial snapshot does not provide a temporal derivative. |
| Pathology and imaging | Register H&E/IHC/IF to coordinates; use blinded annotations for tumor, stroma, necrosis, vessels and immune infiltrates. | Registered pathology and inter-rater agreement. | Game states predict morphology independently of expression-derived proliferation. | Registration artifacts may create spatial associations. |

## Priority 2: complementary diagnostics

| Area | Appropriate use | Falsifiable output | Limitation |
|---|---|---|---|
| Information theory | Conditional mutual information to test whether interactions add information beyond abundance, distance, region and batch. | Incremental information survives conditioning and permutation controls. | Dependence is not direction or causality. |
| Statistical physics | Agent/lattice or field models for clustering, exclusion, interfaces and phase separation; calibrate to spatial summary statistics. | Neutral crowding fails while a specified interaction model reproduces patterns. | Effective phases are not molecular mechanisms without scaling and replication. |
| Machine learning | Compare elastic net, tree models, graph models and interpretable spatial baselines using patient/platform-held-out data. | Pre-specified improvement, calibration and stable attribution. | Predictive importance is not a payoff or causal mechanism. |
| Clinical statistics | Build patient-level features and use Cox/flexible survival models with nested validation and clinical covariates. | External improvement in C-index, time-dependent AUC, Brier score and decision utility. | Prognostic association is not treatment benefit. |
| Decision/control methods from finance | Transfer only generic robust optimization, Bayesian decision theory and stochastic control. Re-derive all biological assumptions. | Robust policy under parameter uncertainty in a biological experiment. | Tumors are not markets; financial analogies are not evidence. |

## Mathematical applicability audit

The merged copy now includes
[`MATHEMATICAL_EXTENSION_REPORT.md`](./MATHEMATICAL_EXTENSION_REPORT.md) and
the reproducible script
[`src/mathematical_applicability_audit.py`](./src/mathematical_applicability_audit.py).
It applies a k-nearest-neighbour graph Laplacian and a truncated spectral zeta
as descriptive diagnostics. The latter is not the Riemann zeta function and
does not estimate a payoff or biological mechanism. The melanoma graph is
disconnected, so its spectrum must not be compared as if it were a single
connected tissue graph.

The audit confirms that graph spectra can describe spatial organization, but
does not change the patient-level conclusion: current game-model superiority
is small, marker-sensitive and uncertain. A Riemann-zeta, Z-transform,
Koopman, stochastic-control or genetic-selection analysis requires a
well-defined observable and, for the dynamic/genetic methods, new data.

The related vorticity proposal is assessed in
[`VORTICITY_APPLICABILITY_REPORT.md`](./VORTICITY_APPLICABILITY_REPORT.md).
The current data have no velocity field or time series, so physical
vorticity `curl(u)` is not identifiable. A curl of marker gradients would be
an arbitrary pseudo-field, not Navier–Stokes vorticity.

The publication-level novelty boundary is recorded in
[`CONTRIBUTION_AND_NOVELTY_SCOPE.md`](./CONTRIBUTION_AND_NOVELTY_SCOPE.md).
The proposed contribution is the auditable separation and validation of
composition, neighbourhood and bilinear interaction terms, not the invention
of spatial graphs, multiscale analysis or evolutionary game theory.

## What should count as evidence for a spatial evolutionary game?

A strong claim requires all of the following:

1. incremental explanatory value over composition and spatial baselines;
2. patient-held-out and spatially robust prediction;
3. stability across deconvolution methods and marker definitions;
4. control for stroma, hypoxia, vessels, immune infiltration, necrosis and
   technical batch;
5. uncertainty propagation into payoff matrices and equilibria;
6. temporal or perturbational prediction;
7. genomic or lineage support for the word "clone", or explicit use of
   "transcriptional state" instead.

The existing in-silico recovery, serial-section concordance and unit tests
support algebraic/software reproducibility. They do not alone establish
causal payoffs, temporal replicator dynamics or clinical treatment benefit.

## Theories deliberately not imported as direct methods

### Quantum mechanics

Molecular biology is quantum at its physical foundation, but the observed
data are aggregate counts, phenotypes and tissue coordinates at scales where
classical stochastic, spatial and statistical models are appropriate. Quantum
games or quantum probability would be metaphors unless a specifically measured
quantum phenomenon were part of the experiment.

### Relativity

Cell migration, diffusion and drug exposure occur at ordinary laboratory
velocities and scales. Relativistic corrections cannot identify spatial
transcriptomic interactions and would add no testable prediction here.

### Millennium Prize problems

These are major pure-mathematics problems, not oncology estimators. Their
solution would not identify payoff parameters or validate a biological
mechanism. Relevant mathematical tools should be selected by the data and
estimand, not by prestige.

### A literal theory of everything

The scientifically coherent objective is a measurement-error-aware,
spatially hierarchical and uncertainty-propagating model with independent
validation. A single universal theory is neither necessary nor identifiable
from this dataset.

## Recommended execution order

1. Validate measurement and deconvolution.
2. Add patient/section hierarchy and spatially appropriate null models.
3. Propagate uncertainty and test marker/graph sensitivity.
4. Add pathology and genomic/lineage evidence.
5. Add longitudinal or perturbation data for dynamics and causal claims.
6. Add survival and adaptive-control analyses only after the previous stages.

## References

- [cell2location](https://doi.org/10.1038/s41587-021-01139-4)
- [RCTD](https://doi.org/10.1038/s41587-021-01148-3)
- [Stereoscope](https://doi.org/10.1038/s41592-020-01047-0)
- [Squidpy](https://doi.org/10.1038/s41592-021-01255-0)
- [Hernán and Robins, Causal Inference](https://www.hsph.harvard.edu/miguel-hernan/causal-inference-book/)
- [Hudgens and Halloran, interference](https://doi.org/10.1093/biomet/95.3.567)
- [Nowak, Evolutionary Dynamics](https://doi.org/10.1007/978-0-691-11774-0)
- [Gatenby et al., adaptive therapy](https://doi.org/10.1371/journal.pone.0010129)
- [Cox proportional hazards model](https://doi.org/10.1111/j.2517-6161.1972.tb00899.x)
- [TRIPOD reporting guideline](https://doi.org/10.1136/bmj.g7594)
