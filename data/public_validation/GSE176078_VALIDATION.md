# GSE176078 patient-held-out validation

The validation uses the six public breast spatial samples audited in
[`GSE176078_AUDIT.md`](./GSE176078_AUDIT.md). The implementation is
[`src/validate_gse176078.py`](../../src/validate_gse176078.py), with numeric
output in [`GSE176078_validation.json`](./GSE176078_validation.json).

## Frozen analysis

- Independent unit: one patient/sample.
- Validation: leave one patient out; train on the other five.
- Spatial radius: 150 µm.
- Baseline: focal marker composition.
- Extended model: focal composition plus mean neighbour composition.
- Outcome: transcriptional proliferation proxy.
- Environment adjustment: hypoxia and vascular modules, with the adjustment
  model fitted on training patients only.
- Excluded spot labels: `Artefact` and `Necrosis`.
- Input: public processed count matrices, not controlled-access FASTQ.

## Result

| Metric | Value |
|---|---:|
| Patients | 6 |
| Mean patient-held-out context-minus-focal R² | -0.06454 |
| Patients with positive increment | 3/6 |
| Bootstrap 95% CI | [-0.39027, +0.18983] |

The interval is wide and includes zero. The result therefore **does not
validate a generalizable improvement** from adding the neighbourhood in this
independent breast cohort. It is a useful negative/uncertain replication
result and should be reported alongside the positive pooled spatial-block
audit rather than omitted.

## Interpretation boundary

The negative result can reflect genuine absence of a transferable signal,
cohort/platform differences, marker mismatch, residual confounding, or low
patient count. It does not prove that spatial context is biologically
irrelevant. Conversely, it cannot support causal communication, evolutionary
payoffs, clonal states, clinical utility, miRNA, or methylome claims.
