# GSE176078 pre-analysis audit

The public processed package was downloaded from the [Zenodo
record](https://zenodo.org/records/4739739) associated with
[GSE176078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078).
The audit is implemented in [`src/audit_gse176078.py`](../../src/audit_gse176078.py)
and its machine-readable output is
[`GSE176078_audit.json`](./GSE176078_audit.json).

## Verified

- Six spatial samples were recovered.
- Each sample maps to one distinct `patientid` in the public metadata:
  `1142243F`, `1160920F`, `CID4290`, `CID4465`, `CID44971`, and `CID4535`.
- Public files include filtered count matrices, spatial coordinates, images and
  spot-level metadata.
- The matrices contain 19,237 or 28,402 genes and 1,127–4,895 filtered spots,
  depending on the sample.
- The required exploratory marker panel is present in every sample:
  tumour-state markers, outcome markers, hypoxia markers and vascular markers.
- Metadata barcodes overlap the public spatial-coordinate barcodes for all
  filtered spots in each sample.

## Access boundary

The GEO record states that raw human sequencing data are deposited under
controlled access. This analysis therefore uses only the public processed
matrices and metadata. No controlled-access account or Google credentials are
needed for this step.

## Still required before claiming validation

- Reconstruct the exact count/coordinate orientation in the model loader.
- Apply consistent spot filtering using the public classification metadata.
- Freeze the marker specification before measuring performance.
- Run patient-held-out validation with all six patients grouped independently.
- Report patient-level effect sizes and uncertainty, not spot-level significance.

The presence of a marker in a matrix does not establish a biological state,
clonal identity, causal interaction, or clinical utility.
