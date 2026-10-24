# Multiscale spatial audit

The existing scan was converted into an explicit audit with
[`src/summarize_multiscale_audit.py`](./src/summarize_multiscale_audit.py).
It evaluates the same tissue at 12 radii, applies Benjamini-Hochberg
correction to exploratory in-sample p-values, and reports spatial-block
out-of-fold (OOF) performance separately.

The output is
[`results/MULTISCALE_AUDIT_SUMMARY.json`](./results/MULTISCALE_AUDIT_SUMMARY.json).

## Results

| Cohort | Correlation length | OOF effect at contact | Best scanned OOF effect | Positive OOF scales |
|---|---:|---:|---:|---:|
| Lung | 475 µm | -0.175 percentage points | +2.680 at 250 µm | 10/12 |
| Breast | 575 µm | +5.462 percentage points | +8.483 at 150 µm | 12/12 |

The lung result is informative: the game model is not better at the contact
scale, but the same-section scan finds a positive effect at intermediate
scales. This is compatible with an interface-scale association, but it is
not an independent validation because the radius was selected after inspecting
the same tissue.

The breast result is also not a causal proof. Its positive OOF result comes
from one original cohort, not eight independent patients. The public
patient-held-out audit remains the stronger generalization check and its
bootstrap interval still includes zero.

## Correct interpretation

- FDR-adjusted in-sample p-values indicate that the radius scan is not merely
  an uncorrected single-test artifact.
- They do **not** establish an evolutionary game mechanism.
- OOF improvement is more relevant than in-sample F-tests, but spatial blocks
  are not independent patients.
- The best radius is exploratory and must be pre-specified and tested in new
  patients or sections.

This is the most suitable additional analysis currently available without
inventing velocity fields or importing unrelated physics. The next decisive
experiment would test the pre-specified 150–300 µm scale range in independent
patients, while preserving the composition-plus-neighborhood baseline.
