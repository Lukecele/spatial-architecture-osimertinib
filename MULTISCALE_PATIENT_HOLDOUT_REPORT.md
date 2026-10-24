# Pre-specified multiscale patient-held-out analysis

The previous same-section scan suggested an interaction scale between
150 and 300 µm. This analysis pre-specifies that window and tests every
radius in it on the 16 paired samples from eight patients in E-MTAB-13530.
For each held-out patient, all remaining patients train both:

- composition plus radius-based neighborhood baseline;
- game interaction model using the same radius-based neighborhood.

The implementation is
[`src/extension_multiscale_patient_holdout.py`](./src/extension_multiscale_patient_holdout.py).
The output is
[`data/extension_pilot/multiscale_patient_holdout_results.json`](./data/extension_pilot/multiscale_patient_holdout_results.json).

## Results

| Radius | Mean patient game-minus-baseline R² | Patients positive | Bootstrap 95% CI |
|---:|---:|---:|---:|
| 150 µm | +0.00185 | 5/8 | [-0.00887, +0.01053] |
| 200 µm | +0.00328 | 5/8 | [-0.00851, +0.01286] |
| 250 µm | +0.00361 | 5/8 | [-0.00811, +0.01331] |
| 300 µm | +0.00531 | 6/8 | [-0.00522, +0.01446] |

All four mean effects are positive, but every bootstrap interval includes
zero. The largest point estimate is at 300 µm, but that does not make 300 µm
a validated optimum. The JSON output is authoritative for the numerical
values. The radius was fixed before this patient-held-out analysis, so no
radius is selected after seeing these results.

## Interpretation rule

Evidence for a robust scale-dependent game contribution would require a
positive mean effect with an interval excluding zero, or a pre-specified
confirmatory test replicated across an independent cohort. A positive result
with an interval covering zero remains exploratory. A negative result would
be evidence against the proposed incremental game contribution at that scale.

This analysis still does not establish causal payoffs, genetic clones,
temporal dynamics, or treatment benefit.
