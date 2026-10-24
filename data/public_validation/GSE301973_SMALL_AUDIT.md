# GSE301973 small-cohort audit

The four-patient EGFR-mutant NSCLC cohort is now locally available from GEO.
Its public biological metadata identify two long-PFS and two short-PFS
pretreatment patients:

- long PFS: 69 and 22 months;
- short PFS: 15 and 12 months.

The public spatial package contains two Visium HD slide objects, with
StarDist-integrated cell-level objects and 2 µm matrices. The audit is
implemented in [`src/audit_gse301973_small.py`](../../src/audit_gse301973_small.py)
and reported numerically in [`GSE301973_small_audit.json`](./GSE301973_small_audit.json).

## What is verified

- the archive contains slide A and slide C;
- cell-level integrated objects contain spatial coordinates and cell-type
  annotations;
- the objects contain `Tumor cell`, airway/alveolar epithelial, immune,
  stromal and endothelial labels;
- EGFR, bypass, EMT, SPP1/CD44, proliferation and stromal markers are present;
- the data are pretreatment, not post-osimertinib.

## Small-study gate

This cohort is suitable for a focused exploratory analysis of spatial programs
associated with pretreatment PFS groups. It is not large enough for a reliable
clinical prediction model or a definitive resistance claim. We will not expand
to high-capacity AI or broad multi-cohort claims unless the small audit first
passes:

1. file and coordinate integrity;
2. patient/slide mapping;
3. cell annotation and QC checks;
4. effect direction consistency under simple models;
5. leave-one-patient-out sensitivity.
