# Osimertinib relevance boundary for GSE189487

The original project focused on lung tumour populations potentially related to
resistance to osimertinib. The public [GSE189487](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189487)
cohort is useful for lung spatial biology, but it is **not** an osimertinib
response cohort:

- it is an untreated LUAD progression atlas;
- it contains spatial samples across AIS, MIA and invasive carcinoma;
- treatment response, drug exposure, progression-free survival and resistance
  labels are not measured in the public spatial matrices.

The audit in [`src/lung_osimertinib_relevance_audit.py`](../../src/lung_osimertinib_relevance_audit.py)
therefore evaluates only the spatial organisation of candidate transcriptional
programs:

- EGFR/bypass signalling;
- MET/AXL/FGFR/IGF1R-associated bypass;
- EMT-associated state;
- stem-like markers;
- proliferation;
- epithelial/alveolar identity.

The numeric output is
[`GSE189487_OSIMERTINIB_RELEVANCE.json`](./GSE189487_OSIMERTINIB_RELEVANCE.json).
Moran-style spatial summaries are descriptive and are not drug-response
statistics.

## What would be required for a true osimertinib validation

- EGFR-mutant tumours with documented osimertinib exposure;
- pre-treatment and on-treatment or progression samples;
- response/progression labels;
- patient-level grouping;
- ideally matched DNA/CNV and single-cell or spatial data;
- independent patient-held-out validation.

Until then, terms such as “resistant population” must be written as
“candidate resistance-associated transcriptional state”.
