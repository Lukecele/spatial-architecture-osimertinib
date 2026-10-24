# GSE189487 lung-adapted validation

The analysis uses the public processed matrices from
[GSE189487](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189487).
The implementation is [`src/validate_gse189487.py`](../../src/validate_gse189487.py)
and the machine-readable output is
[`GSE189487_validation.json`](./GSE189487_validation.json).

## Universal skeleton

- focal composition versus focal composition plus spatial neighbourhood;
- 150 µm radius;
- leave-one-spatial-sample-out validation;
- patient/sample-level effect reporting;
- environmental adjustment fitted on training samples only.

## Lung-specific adaptation

The marker representation was not copied from breast. It used four
lung-relevant descriptive compartments:

- epithelial/alveolar: `EPCAM`, `KRT8`, `SCGB1A1`, `SFTPC`;
- invasive/mesenchymal: `VIM`, `AXL`, `CDH2`, `KRT19`;
- CAF: `COL1A1`, `COL3A1`, `ACTA2`, `FAP`;
- immune: `CD68`, `CD3D`, `CD8A`, `MS4A1`, `NKG7`.

The outcome proxy was `MKI67`, `PCNA`, `TOP2A`, `CDK1`; hypoxia and vascular
modules were used as environmental covariates.

## Result

| Metric | Value |
|---|---:|
| Spatial samples | 6 |
| Mean held-out context-minus-focal R² | -0.04793 |
| Samples with positive increment | 4/6 |
| Bootstrap 95% CI | [-0.29986, +0.14579] |

This is a negative/uncertain replication of a general predictive improvement.
The interval includes zero and is wide. It does not invalidate the universal
skeleton by itself: the cohort is small, stage and sample are not fully
separable, and this is still a transcriptomic proxy analysis. It does show that
lung-specific adaptation alone did not produce a robust held-out gain at the
pre-specified 150 µm radius.

## Boundary

The six public spatial samples are `TD1`, `TD2`, `TD3`, `TD5`, `TD6`, and `TD8`.
The GEO study describes nine patients overall and six paired spatial samples;
the exact stage/sample mapping must be retained from the study metadata before
any stage-specific claim. No causal interaction, clonal identity, payoff,
clinical utility, miRNA or methylome conclusion follows from this result.
