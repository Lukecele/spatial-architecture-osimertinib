# InSituGame extended pilot

This directory is an independent copy of the original InSituGame project.
The original project at `/home/luca/spatial-game-assay` was not modified.

## Scope

The extension integrates a small, paired technical pilot from the public
E-MTAB-13530 NSCLC study:

- `P10_B1`: adjacent tissue from Patient 10;
- `P10_T1`: tumour tissue from Patient 10.

Both samples were processed with the existing InSituGame implementation. No
changes were made to the assay, marker definitions, estimator, Nash solver,
or validation code.

## Integrated artifacts

The input and output files are under
[`data/extension_pilot/`](./data/extension_pilot/):

- `pilot_P10_B1/` and `pilot_P10_T1/`: filtered HDF5 matrices and spatial
  files;
- [`metadata_audit.tsv`](./data/extension_pilot/metadata_audit.tsv): cohort
  metadata audit;
- [`pilot_manifest.tsv`](./data/extension_pilot/pilot_manifest.tsv): patient,
  tissue and matrix dimensions;
- [`pilot_qc.tsv`](./data/extension_pilot/pilot_qc.tsv): minimal QC summary;
- [`pilot_assay_results.json`](./data/extension_pilot/pilot_assay_results.json):
  exploratory assay output.

## Results

| Sample | Tissue | InSituGame regime | In-sample R² |
|---|---|---|---:|
| `P10_B1` | adjacent | Resistant Dominance | 0.0151 |
| `P10_T1` | tumour | Stag Hunt (Bistable), `x_sep = 0.206` | 0.0201 |

These are exploratory paired results. They do not establish a biological
regime, clinical generalization, causal payoff interpretation, or superiority
of the method. The low R² values and the different library depths require
additional QC and replication before biological interpretation.

## Baseline comparison

The script [`src/extension_pilot_baselines.py`](./src/extension_pilot_baselines.py)
compares the existing game design with a mean-only baseline and a
composition-plus-neighborhood baseline. The comparison is in-sample and is
not a test of external validity.

| Sample | Composition + neighborhood R² | Game R² | Increment |
|---|---:|---:|---:|
| `P10_B1` | 0.0139 | 0.0151 | 0.0012 |
| `P10_T1` | 0.0168 | 0.0201 | 0.0034 |

The game design adds only a small in-sample increase in this two-sample pilot.
This is useful as a falsifiable baseline check, but it does not establish that
the interaction terms are biologically necessary. Patient-held-out validation
and spatially structured nulls remain required.

The machine-readable output is
[`pilot_baseline_comparison.json`](./data/extension_pilot/pilot_baseline_comparison.json).

## Spatial block validation and null model

The same script also performs five-fold spatial block cross-validation and a
100-permutation shuffle of neighborhood frequencies. Results:

| Sample | Composition block-CV R² | Game block-CV R² | Shuffle-null empirical p |
|---|---:|---:|---:|
| `P10_B1` | 0.0089 | 0.0061 | 0.0198 |
| `P10_T1` | -0.0048 | -0.0101 | 0.0099 |

The negative held-out R² values in `P10_T1` and the lower game score than the
composition baseline mean that this pilot does **not** demonstrate spatial
generalization. The shuffle p-values are not sufficient to rescue that
conclusion: the null preserves neither all spatial dependence nor patient
heterogeneity, and the observed score is in-sample. They are diagnostic only.

## Leave-one-patient-out extension

To move beyond the single-patient pilot, filtered matrices and spatial files
were added for paired samples from Patients 11 and 15. A leave-one-patient-out
audit now trains on the other two patients and tests on the held-out patient,
with composition-plus-neighborhood regression as the baseline.

| Held-out sample | Baseline R² | Game R² | Game minus baseline |
|---|---:|---:|---:|
| `P10_B1` | -0.0149 | 0.0031 | +0.0180 |
| `P10_T1` | -0.0010 | -0.0202 | -0.0191 |
| `P11_B1` | 0.0023 | 0.0059 | +0.0037 |
| `P11_T1` | 0.0293 | 0.0905 | +0.0612 |
| `P15_B1` | -0.0186 | -0.0134 | +0.0052 |
| `P15_T1` | -0.0173 | -0.0069 | +0.0105 |

This is still a small, same-study exploratory cohort, not clinical validation.
The game model is not uniformly superior and several held-out R² values are
negative. The result supports continued patient-level testing while ruling
out a claim of universal generalization.

Machine-readable output:
[`patient_holdout_results.json`](./data/extension_pilot/patient_holdout_results.json).

The analysis was then extended to all eight patients with paired files
available in this study record: Patients 10, 11, 15, 16, 17, 19, 24 and 25.
Patient-level mean game-minus-baseline R² was `+0.00185`; five of eight
patients had a positive mean difference. A 10,000-replicate patient bootstrap
gave a 95% interval of `[-0.00883, +0.01050]`, which includes zero.

This is the strongest current result from the extension and it is still
negative for a claim of reliable superiority: the average advantage is small,
uncertain and heterogeneous. The patient-level summary is in
[`patient_holdout_summary.json`](./data/extension_pilot/patient_holdout_summary.json).

## Marker sensitivity

Four marker panels were evaluated across all 16 samples using the same
five-fold spatial block procedure:

| Marker panel | Mean baseline block-CV R² | Mean game block-CV R² | Samples where game wins |
|---|---:|---:|---:|
| Reference | 0.0275 | 0.0249 | 4/16 |
| Remove VIM from T_R | 0.0310 | 0.0277 | 5/16 |
| Expanded epithelial T_S | 0.0300 | 0.0284 | 7/16 |
| Expanded CAF | 0.0263 | 0.0237 | 4/16 |

The game model is below the composition-plus-neighborhood baseline on mean
held-out R² for all four panels. This is a negative robustness result: the
current evidence does not show that the apparent game signal is stable or
incremental under marker redefinition.

Output:
[`marker_sensitivity_results.json`](./data/extension_pilot/marker_sensitivity_results.json).

## Sequencing-depth confounding audit

The QC table was rebuilt for all 16 samples and joined to the patient-held-out
results. The correlation between the game-minus-baseline R² difference and
sample median depth was:

- median UMI: `r = -0.426`;
- median detected features: `r = -0.446`.

These are diagnostic correlations, not adjusted effects. They indicate that
technical depth may be entangled with the apparent model advantage and must be
included in a future hierarchical measurement-error model. They do not prove
that depth causes the result, because tissue type and patient are also
structured with the measurements.

Output:
[`depth_confounding_audit.json`](./data/extension_pilot/depth_confounding_audit.json).

## Hierarchical QC-adjusted audit

Because a full random-effects package is not installed in the current
environment, a dependency-free audit was implemented instead. It preserves
patient-held-out evaluation, adds log median UMI and log detected features as
sample-level QC covariates, and summarizes performance at the patient level.
It is explicitly not a substitute for a Bayesian or mixed-effects model.

The QC-adjusted mean patient-level game-minus-baseline R² was `+0.00395`, with
a 10,000-replicate patient bootstrap 95% interval of
`[-0.01289, +0.01994]`. The interval still includes zero and is wider than
the unadjusted estimate. This does not rescue the incremental game claim; it
shows that QC adjustment must be part of the final model and uncertainty
reporting.

Output:
[`hierarchical_qc_audit.json`](./data/extension_pilot/hierarchical_qc_audit.json).

## Validation status

The original project test suite was run with:

```bash
python3 -m unittest discover -s tests -v
```

Result: 19 tests passed.

The pilot is therefore integrated as a reproducible evidence layer, not as a
replacement for the original results. Public literature datasets may be used
as further technical fixtures or independent transfer analyses, but they must
remain separated by patient, tissue, platform and batch.
