# Core validated model

## Scope

The defensible core of the merged study is a spatial-context assay for
spatial transcriptomics. It is not, at this stage, a validated tumour
evolutionary game or a treatment-control model.

## Scientific question

> Does measured neighbourhood context add reproducible information about a
> molecular outcome beyond the composition measured at the focal spot?

## Core inputs

- spot coordinates;
- gene-level counts or normalized expression;
- a pre-specified marker panel;
- a molecular outcome proxy;
- patient, section and batch identifiers when available.

## Core transformations

1. Normalize expression with a documented library-size rule.
2. Compute pre-specified module scores.
3. Represent spot composition on a simplex only as a descriptive state
   representation.
4. Construct a documented neighbourhood graph or fixed-radius field.
5. Compute focal composition and neighbourhood composition.
6. Fit transparent baseline models before any interaction model.
7. Evaluate spatial-block and patient-held-out performance.
8. Aggregate uncertainty at the biological replicate level.

## Model hierarchy

| Model | Role |
|---|---|
| Outcome mean | Minimal reference |
| Focal composition | Tests spot-level information |
| Focal plus neighbourhood composition | Tests spatial context |
| Bilinear interaction | Candidate exploratory extension |

The bilinear model is not automatically called a game model in the primary
analysis. Its coefficients are conditional regression parameters unless
independent biological and temporal evidence supports a stronger
interpretation.

## Acceptance criteria

An interaction claim can be promoted beyond exploratory status only if:

1. the interaction model improves a pre-specified spatial baseline;
2. the improvement survives held-out spatial blocks;
3. the result is evaluated at patient/section level rather than spot level;
4. uncertainty excludes a practically negligible effect;
5. marker and graph choices do not reverse the conclusion;
6. technical depth and batch do not explain the gain;
7. the result replicates in an independent cohort or perturbation.

Failure at any of these gates means the interaction remains exploratory.

## Current status

The software and mathematical invariants are verified. Spatial structure is
present. The patient-held-out multiscale effects from 150--300 µm are small
and all bootstrap intervals include zero. Therefore:

- the spatial-context assay is a valid reproducible analysis framework;
- the bilinear/game extension is not yet a validated biological effect;
- the published microRNA/methylome mechanisms are not measurable in the
  current cache;
- causal, genetic and temporal claims remain outside the core.

## Modules kept outside the core

- Nash classification and replicator dynamics;
- adaptive therapy and Pontryagin control;
- miRNA and methylome mechanisms without matched measurements;
- vorticity, quantum and relativistic analogies;
- clinical decision claims.

These modules may remain as hypotheses or simulations, but must not define the
primary conclusion of the study.

## Reproducibility rule

Any future extension must report, in order:

1. data availability and provenance;
2. pre-specified biological question;
3. baseline;
4. held-out split;
5. uncertainty unit;
6. sensitivity results;
7. interpretation boundary.

This prevents model complexity from being mistaken for validation.

The executable minimal audit is
[`src/run_core_spatial_context_audit.py`](./src/run_core_spatial_context_audit.py),
with machine-readable output in
[`results/CORE_SPATIAL_CONTEXT_AUDIT.json`](./results/CORE_SPATIAL_CONTEXT_AUDIT.json).

The public patient-held-out implementation is
[`src/core_patient_holdout.py`](./src/core_patient_holdout.py), with results
in [`data/extension_pilot/core_patient_holdout.json`](./data/extension_pilot/core_patient_holdout.json).
