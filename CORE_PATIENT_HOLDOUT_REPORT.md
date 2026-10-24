# Core patient-held-out validation

This analysis applies the minimal non-game model to the eight paired public
patients. The radius is fixed at 150 µm, matching the core spatial audit.
For each held-out patient, models are trained on the other seven patients:

- focal composition only;
- focal composition plus neighbourhood composition.

The implementation is
[`src/core_patient_holdout.py`](./src/core_patient_holdout.py), and the
machine-readable output is
[`data/extension_pilot/core_patient_holdout.json`](./data/extension_pilot/core_patient_holdout.json).

## Results

| Metric | Result |
|---|---:|
| Patients | 8 |
| Mean patient context-minus-focal R² | +0.00091 |
| Patients with positive effect | 6/8 |
| Bootstrap 95% CI | [-0.00163, +0.00316] |

## Boundary

The mean improvement is small and its patient-level bootstrap interval includes
zero. Thus this analysis does not establish a reproducible patient-level
increment beyond focal composition. It is still useful as a transparent
negative/uncertain validation result: the large spatial-block improvements
observed in the pooled core audit do not translate into a clearly separated
patient-held-out effect in this eight-patient public extension.

The cohort is one public study, the outcome is a molecular proxy, and the input
is restricted to a marker panel. It does not validate causal communication,
genetic clones, miRNA, methylome, therapy response, or clinical utility.
