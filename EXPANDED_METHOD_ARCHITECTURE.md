# Expanded method architecture

## Strategic change

The project should not require an evolutionary-game interpretation to succeed.
The transferable contribution is a reproducible way to test whether local
spatial context adds information beyond the focal spot, with tumour-specific
biological representations.

The game-theoretic model is therefore an optional hypothesis module, not the
scientific foundation.

## Universal skeleton

For each tumour and dataset:

1. define a biologically justified representation of the focal spot;
2. define an outcome measured independently from that representation;
3. construct one or more spatial context operators;
4. compare focal-only, context-augmented, and null models;
5. validate at the patient/donor level;
6. quantify uncertainty and sensitivity;
7. interpret only the evidence supported by the data.

The skeleton is reusable. Marker panels, spatial scales, graph construction,
outcomes, and nuisance covariates remain tumour- and platform-specific.

## Model families to compare

### A. Spatial statistics

- Moran's I and Geary's C for spatial autocorrelation;
- local indicators of spatial association;
- variograms and distance-decay curves;
- permutation or restricted-null tests.

Purpose: establish whether the signal is spatial before adding biological
mechanisms.

### B. Graph and signal models

- radius and k-nearest-neighbour graphs;
- graph Laplacian smoothing as a diagnostic;
- graph signal energy and local heterogeneity;
- graph-regularized regression.

Purpose: quantify geometry and neighbourhood structure without calling graph
coefficients biological payoffs.

### C. Multiscale context

- several pre-specified radii;
- distance-weighted neighbourhoods;
- scale-response curves;
- patient-held-out selection or externally fixed scale.

Purpose: determine whether the relevant context is local, intermediate, or
regional.

### D. Cell-composition and communication models

- deconvolution or cell-type abundance where reference data exist;
- ligand-receptor analysis as an exploratory association layer;
- composition-aware spatial regression;
- explicit separation of co-localisation from communication.

Purpose: move from generic expression modules toward interpretable cell
neighbourhoods without claiming causality.

### E. Hierarchical and predictive models

- patient/section random effects;
- nested cross-validation;
- grouped permutation tests;
- calibrated uncertainty;
- benchmark against simple smoothing and non-spatial baselines.

Purpose: test reproducibility and prevent spot-level pseudoreplication.

### F. Imaging and pathology integration

- histology-informed regions;
- tumour/stroma/necrosis boundaries;
- morphology features;
- pathology-annotated validation.

Purpose: distinguish molecular neighbourhood from anatomical or technical
structure.

### G. Game-theoretic module

Retain only as an optional extension. Its claims require evidence for:

- state definitions;
- interaction terms;
- temporal or perturbational dynamics;
- stable parameters across patients;
- biological interpretation of coefficients.

Without these, report the module as a nonlinear interaction model, not as
validated evolutionary payoff dynamics.

## Explicitly excluded claims

The expanded architecture does not by itself establish:

- causality;
- clonal identity or lineage;
- evolutionary fitness;
- therapeutic response;
- miRNA or methylome mechanisms;
- physical vorticity or fluid dynamics.

Those require additional measurements or interventions.

## Recommended execution order

1. Complete dataset-level audits and patient mapping.
2. Establish spatial autocorrelation and restricted nulls.
3. Compare focal-only, smoothing, graph, and multiscale models.
4. Add tumour-specific cell-composition representations.
5. Add pathology and imaging where available.
6. Evaluate communication hypotheses as exploratory associations.
7. Revisit game-theoretic interactions only if simpler models leave a stable,
   reproducible residual interaction.

The broader inventory and filtering gates are specified in
[`METHOD_INVENTORY_AND_FILTER.md`](./METHOD_INVENTORY_AND_FILTER.md), with a
machine-readable tier registry in [`METHOD_REGISTRY.json`](./METHOD_REGISTRY.json).

## Core claim to target

> A tumour-adapted, multiscale spatial-context protocol can quantify whether
> neighbouring tissue adds reproducible information beyond the focal molecular
> state, while separating geometric association, predictive utility, and
> biological mechanism.
