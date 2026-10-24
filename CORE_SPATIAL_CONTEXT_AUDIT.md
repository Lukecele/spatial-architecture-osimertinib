# Minimal non-game spatial-context audit

The executable audit in
[`src/run_core_spatial_context_audit.py`](./src/run_core_spatial_context_audit.py)
uses the bundled gene-level cache and compares:

1. focal spot composition;
2. focal composition plus mean neighbourhood composition within 150 µm.

The output is
[`results/CORE_SPATIAL_CONTEXT_AUDIT.json`](./results/CORE_SPATIAL_CONTEXT_AUDIT.json).
The outcome is the proliferation module after linear adjustment for hypoxia and
vascular modules. This is a descriptive molecular proxy, not a growth rate.

## Results

| Cohort | Focal OOF R² | Focal + neighbourhood OOF R² | Increment |
|---|---:|---:|---:|
| Breast | 0.1074 | 0.1573 | +0.0498 |
| Melanoma | 0.1150 | 0.1937 | +0.0787 |
| Lung | 0.2236 | 0.2489 | +0.0249 |

At this fixed radius, neighbourhood composition improves spatial-block OOF
prediction over focal composition in all three cached cohorts. This validates
the narrower statement that measured spatial context contains reproducible
predictive information for this molecular proxy within these sections.

The audit also permutes neighbourhood composition across spots 100 times and
recomputes the spatial-block OOF increment. The empirical null is reported per
cohort in the JSON output. This tests whether the observed gain is compatible
with arbitrary reassignment of neighbourhood features; it is not a substitute
for patient-level replication or a causal null.

Observed versus permutation results:

| Cohort | Observed OOF increment | Null mean | Null maximum | Empirical p |
|---|---:|---:|---:|---:|
| Breast | +0.04984 | -0.00013 | +0.00151 | 0.0099 |
| Melanoma | +0.07868 | -0.00035 | +0.00185 | 0.0099 |
| Lung | +0.02531 | -0.00036 | +0.00048 | 0.0099 |

This strengthens the narrow statement that the measured neighbourhood
features carry non-random spatial-context information within these sections.
It still does not turn the result into patient-level, causal, or game-theoretic
validation.

It does **not** validate:

- a bilinear interaction or evolutionary game;
- causal cell-cell communication;
- microRNA or methylome mechanisms;
- genetic clones;
- patient-level generalization, because these cached cohorts are sections and
  spatial blocks are not independent patients.

The result supports making the non-game spatial-context model primary and
keeping game terms as a secondary hypothesis.
